#include "io_service_pool.hpp"

io_service_pool::io_service_pool(const unsigned concurrency) : w(new work(*this))
{
	reserve(concurrency);
	for (int i = 0; i < concurrency; ++i)
	{
		emplace_back(async(launch::async, [&]()
		{
			run();
		}));
	}
}

void io_service_pool::wait()
{
	w.reset();
	for (auto& f : *this)
	{
		f.get();
	}
}
