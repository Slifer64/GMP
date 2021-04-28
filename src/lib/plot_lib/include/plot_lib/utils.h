#ifndef AS64_QTPLOT_UTILS_H
#define AS64_QTPLOT_UTILS_H

#include <mutex>
#include <condition_variable>

namespace as64_
{

namespace pl_
{

class Semaphore
{
public:
  void notify();
  void wait();
  bool try_wait();

private:
  std::mutex mutex_;
  std::condition_variable condition_;
  // unsigned long count_ = 0; // Initialized as locked.
  bool count_ = false;  // Initialized as locked.
};

} // namespace pl_

} // namespace as64_

#endif // AS64_QTPLOT_UTILS_H
