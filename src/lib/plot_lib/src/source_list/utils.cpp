#include <plot_lib/utils.h>


namespace as64_
{

namespace pl_
{

void Semaphore::notify()
{
  std::lock_guard<decltype(mutex_)> lock(mutex_);
  // ++count_;
  count_ = true;
  condition_.notify_one();
}

void Semaphore::wait()
{
  std::unique_lock<decltype(mutex_)> lock(mutex_);

  // Handle spurious wake-ups.
  while(!count_) condition_.wait(lock);
  // --count_;
  count_ = false;
}

bool Semaphore::try_wait()
{
  std::lock_guard<decltype(mutex_)> lock(mutex_);
  if(count_)
  {
    // --count_;
    count_ = false;
    return true;
  }
  return false;
}

} // namespace pl_

} // namespace as64_
