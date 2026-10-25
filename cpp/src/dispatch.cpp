#include "dispatch.hpp"
#include "worker.hpp"
#include <mpi.h>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <mutex>

namespace genomic {

namespace {

using Clock = std::chrono::steady_clock;

struct WorkerState {
  NodeInfo node;
  int64_t jobs = 0;
  int64_t failures = 0;
  int active = 0;
  TaskType current = TaskType::kShutdown;
  double compute_ms = 0;
  int64_t bytes_sent = 0;
  int64_t bytes_received = 0;
};

int g_world_size = 1;
bool g_thread_multiple = false;

// Guards g_workers and g_idle. Never held during MPI calls or computation.
std::mutex g_mu;
std::condition_variable g_idle_cv;
std::deque<int> g_idle;              // ranks waiting for work (FIFO = round robin)
std::vector<WorkerState> g_workers;  // indexed by rank

// Serializes MPI when the library lacks MPI_THREAD_MULTIPLE.
std::mutex g_mpi_serial;

double ms_since(Clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

int acquire_worker(TaskType type) {
  std::unique_lock<std::mutex> lock(g_mu);
  g_idle_cv.wait(lock, [] { return !g_idle.empty(); });
  int w = g_idle.front();
  g_idle.pop_front();
  g_workers[w].active = 1;
  g_workers[w].current = type;
  return w;
}

// Send one framed task to worker w and wait for its framed reply.
std::string mpi_round_trip(int w, const std::string& task) {
  std::unique_lock<std::mutex> serial(g_mpi_serial, std::defer_lock);
  if (!g_thread_multiple) serial.lock();

  MPI_Send(task.data(), static_cast<int>(task.size()), MPI_BYTE,
           w, kTagTask, MPI_COMM_WORLD);

  MPI_Message handle;
  MPI_Status status;
  MPI_Mprobe(w, kTagResult, MPI_COMM_WORLD, &handle, &status);
  int len = 0;
  MPI_Get_count(&status, MPI_BYTE, &len);
  std::string reply(static_cast<size_t>(len), '\0');
  MPI_Mrecv(reply.data(), len, MPI_BYTE, &handle, MPI_STATUS_IGNORE);
  return reply;
}

bool run_local(TaskType type, const std::string& payload, std::string* result,
               DispatchInfo* info, std::string* error, Clock::time_point t0) {
  {
    std::lock_guard<std::mutex> lock(g_mu);
    ++g_workers[0].active;
    g_workers[0].current = type;
  }
  auto c0 = Clock::now();
  bool ok = execute_task(type, payload, result);
  double compute_ms = ms_since(c0);
  {
    std::lock_guard<std::mutex> lock(g_mu);
    auto& s = g_workers[0];
    --s.active;
    ++s.jobs;
    if (!ok) ++s.failures;
    s.compute_ms += compute_ms;
  }
  info->worker_rank = 0;
  info->compute_ms = compute_ms;
  info->total_ms = ms_since(t0);
  if (!ok) *error = *result;
  return ok;
}

}  // namespace

void dispatch_init(int world_size, bool thread_multiple, std::vector<NodeInfo> nodes) {
  std::lock_guard<std::mutex> lock(g_mu);
  g_world_size = world_size;
  g_thread_multiple = thread_multiple;
  g_workers.assign(static_cast<size_t>(world_size), WorkerState{});
  for (int r = 0; r < world_size; ++r) {
    g_workers[r].node = r < static_cast<int>(nodes.size()) ? nodes[r] : NodeInfo{};
    g_workers[r].node.rank = r;
  }
  g_idle.clear();
  for (int r = 1; r < world_size; ++r) g_idle.push_back(r);
}

int dispatch_worker_count() {
  return (g_world_size > 1) ? (g_world_size - 1) : 0;
}

bool dispatch_task(TaskType type, const std::string& payload,
                   std::string* result, DispatchInfo* info, std::string* error) {
  auto t0 = Clock::now();
  if (g_world_size <= 1) return run_local(type, payload, result, info, error, t0);

  int w = acquire_worker(type);
  std::string task = frame_task(type, payload);
  std::string reply = mpi_round_trip(w, task);

  int32_t status = 1;
  int64_t compute_us = 0;
  bool framed = unframe_result(reply, &status, &compute_us, result);
  bool ok = framed && status == 0;
  {
    std::lock_guard<std::mutex> lock(g_mu);
    auto& s = g_workers[w];
    s.active = 0;
    ++s.jobs;
    if (!ok) ++s.failures;
    s.compute_ms += compute_us / 1000.0;
    s.bytes_sent += static_cast<int64_t>(task.size());
    s.bytes_received += static_cast<int64_t>(reply.size());
    g_idle.push_back(w);
  }
  g_idle_cv.notify_one();

  info->worker_rank = w;
  info->compute_ms = compute_us / 1000.0;
  info->total_ms = ms_since(t0);
  if (!ok) *error = framed ? *result : "malformed reply from worker";
  return ok;
}

std::vector<WorkerSnapshot> dispatch_snapshot() {
  std::lock_guard<std::mutex> lock(g_mu);
  std::vector<WorkerSnapshot> out;
  int first = g_world_size > 1 ? 1 : 0;
  for (int r = first; r < static_cast<int>(g_workers.size()); ++r) {
    const auto& s = g_workers[r];
    WorkerSnapshot snap;
    snap.node = s.node;
    snap.jobs = s.jobs;
    snap.failures = s.failures;
    snap.active = s.active;
    snap.current_task = s.active > 0 ? task_name(s.current) : "";
    snap.compute_ms = s.compute_ms;
    snap.bytes_sent = s.bytes_sent;
    snap.bytes_received = s.bytes_received;
    out.push_back(std::move(snap));
  }
  return out;
}

void dispatch_shutdown_workers() {
  std::string msg = frame_task(TaskType::kShutdown, "");
  for (int w = 1; w < g_world_size; ++w) {
    MPI_Send(msg.data(), static_cast<int>(msg.size()), MPI_BYTE,
             w, kTagTask, MPI_COMM_WORLD);
  }
}

}  // namespace genomic
