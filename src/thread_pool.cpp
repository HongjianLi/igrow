#include "thread_pool.hpp"

thread_pool::thread_pool(const size_t num_threads) : num_scheduled_tasks(0), num_completed_tasks(0), exiting(false)
{
	reserve(num_threads);
	for (size_t thread_id = 0; thread_id < num_threads; ++thread_id)
	{
		push_back(thread([&,thread_id]()
		{
			while (true)
			{
				{
					unique_lock<mutex> lock(m);
					if (num_scheduled_tasks == tasks.size())
					{
						task_incoming.wait(lock);
					}
				}
				if (exiting) return;
				while (true)
				{
					size_t task_id;
					{
						unique_lock<mutex> lock(m);
						if (num_scheduled_tasks == tasks.size()) break;
						task_id = num_scheduled_tasks++;
					}
					tasks[task_id].operator()();
					{
						unique_lock<mutex> lock(m);
						++num_completed_tasks;
					}
					task_completion.notify_one();
				}
			}
		}));
	}
}

void thread_pool::enqueue(packaged_task<void()>&& task)
{
	unique_lock<mutex> lock(m);
	tasks.push_back(static_cast<packaged_task<void()>&&>(task));
	task_incoming.notify_one();
}

void thread_pool::synchronize()
{
	unique_lock<mutex> lock(m);
	while (num_completed_tasks < tasks.size())
	{
		task_completion.wait(lock);
	}
	for (auto& task : tasks)
	{
		task.get_future().get();
	}
	tasks.clear();
	num_scheduled_tasks = 0;
	num_completed_tasks = 0;
}

thread_pool::~thread_pool()
{
	exiting = true;
	task_incoming.notify_all();
	for (auto& thread : *this)
	{
		thread.join();
	}
}
