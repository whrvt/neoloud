/*
SoLoud audio engine
Copyright (c) 2013-2015 Jari Komppa

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.
*/

#if __cplusplus >= 201703L
#include <chrono>
#include <mutex>
#include <thread>
#else
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#if defined(__MINGW32__) || defined(__MINGW64__)
#include <pthread.h>
#endif
#else
#include <inttypes.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#endif
#endif

#include "soloud_thread.h"

namespace SoLoud
{
namespace Thread
{

#if __cplusplus >= 201703L

struct ThreadHandleData
{
	std::thread thread;
};

void *createMutex()
{
	return new std::mutex();
}

void destroyMutex(void *aHandle)
{
	std::mutex *mutex = static_cast<std::mutex *>(aHandle);
	delete mutex;
}

void lockMutex(void *aHandle)
{
	std::mutex *mutex = static_cast<std::mutex *>(aHandle);
	if (mutex)
	{
		mutex->lock();
	}
}

void unlockMutex(void *aHandle)
{
	std::mutex *mutex = static_cast<std::mutex *>(aHandle);
	if (mutex)
	{
		mutex->unlock();
	}
}

struct soloud_thread_data
{
	threadFunction mFunc;
	void *mParam;
};

static void threadfunc(soloud_thread_data *d)
{
	d->mFunc(d->mParam);
	delete d;
}

ThreadHandle createThread(threadFunction aThreadFunction, void *aParameter)
{
	soloud_thread_data *d = new soloud_thread_data;
	d->mFunc = aThreadFunction;
	d->mParam = aParameter;

	ThreadHandleData *threadHandle = new ThreadHandleData;
	threadHandle->thread = std::thread(threadfunc, d);
	return threadHandle;
}

void sleep(int aMSec)
{
	std::this_thread::sleep_for(std::chrono::milliseconds(aMSec));
}

void wait(ThreadHandle aThreadHandle)
{
	if (aThreadHandle->thread.joinable())
	{
		aThreadHandle->thread.join();
	}
}

void release(ThreadHandle aThreadHandle)
{
	delete aThreadHandle;
}

int getTimeMillis()
{
	auto now = std::chrono::steady_clock::now();
	auto duration = now.time_since_epoch();
	auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
	return static_cast<int>(millis);
}

#else

#ifdef WINDOWS_VERSION
struct ThreadHandleData
{
	HANDLE thread;
};

void *createMutex()
{
	CRITICAL_SECTION *cs = new CRITICAL_SECTION;
	InitializeCriticalSectionAndSpinCount(cs, 100);
	return (void *)cs;
}

void destroyMutex(void *aHandle)
{
	CRITICAL_SECTION *cs = (CRITICAL_SECTION *)aHandle;
	DeleteCriticalSection(cs);
	delete cs;
}

void lockMutex(void *aHandle)
{
	CRITICAL_SECTION *cs = (CRITICAL_SECTION *)aHandle;
	if (cs)
	{
		EnterCriticalSection(cs);
	}
}

void unlockMutex(void *aHandle)
{
	CRITICAL_SECTION *cs = (CRITICAL_SECTION *)aHandle;
	if (cs)
	{
		LeaveCriticalSection(cs);
	}
}

struct soloud_thread_data
{
	threadFunction mFunc;
	void *mParam;
};

static DWORD WINAPI threadfunc(LPVOID d)
{
	soloud_thread_data *p = (soloud_thread_data *)d;
	p->mFunc(p->mParam);
	delete p;
	return 0;
}

ThreadHandle createThread(threadFunction aThreadFunction, void *aParameter)
{
	soloud_thread_data *d = new soloud_thread_data;
	d->mFunc = aThreadFunction;
	d->mParam = aParameter;
	HANDLE h = CreateThread(NULL, 0, threadfunc, d, 0, NULL);
	if (0 == h)
	{
		return 0;
	}
	ThreadHandleData *threadHandle = new ThreadHandleData;
	threadHandle->thread = h;
	return threadHandle;
}

void sleep(int aMSec)
{
	Sleep(aMSec);
}

void wait(ThreadHandle aThreadHandle)
{
	WaitForSingleObject(aThreadHandle->thread, INFINITE);
}

void release(ThreadHandle aThreadHandle)
{
	CloseHandle(aThreadHandle->thread);
	delete aThreadHandle;
}

int getTimeMillis()
{
	return GetTickCount();
}

#else // pthreads
struct ThreadHandleData
{
	pthread_t thread;
};

void *createMutex()
{
	pthread_mutex_t *mutex;
	mutex = new pthread_mutex_t;

	pthread_mutexattr_t attr;
	pthread_mutexattr_init(&attr);

	pthread_mutex_init(mutex, &attr);

	return (void *)mutex;
}

void destroyMutex(void *aHandle)
{
	pthread_mutex_t *mutex = (pthread_mutex_t *)aHandle;

	if (mutex)
	{
		pthread_mutex_destroy(mutex);
		delete mutex;
	}
}

void lockMutex(void *aHandle)
{
	pthread_mutex_t *mutex = (pthread_mutex_t *)aHandle;
	if (mutex)
	{
		pthread_mutex_lock(mutex);
	}
}

void unlockMutex(void *aHandle)
{
	pthread_mutex_t *mutex = (pthread_mutex_t *)aHandle;
	if (mutex)
	{
		pthread_mutex_unlock(mutex);
	}
}

struct soloud_thread_data
{
	threadFunction mFunc;
	void *mParam;
};

static void *threadfunc(void *d)
{
	soloud_thread_data *p = (soloud_thread_data *)d;
	p->mFunc(p->mParam);
	delete p;
	return 0;
}

ThreadHandle createThread(threadFunction aThreadFunction, void *aParameter)
{
	soloud_thread_data *d = new soloud_thread_data;
	d->mFunc = aThreadFunction;
	d->mParam = aParameter;

	ThreadHandleData *threadHandle = new ThreadHandleData;
	pthread_create(&threadHandle->thread, NULL, threadfunc, (void *)d);
	return threadHandle;
}

void sleep(int aMSec)
{
	// usleep(aMSec * 1000);
	struct timespec req = {0};
	req.tv_sec = 0;
	req.tv_nsec = aMSec * 1000000L;
	nanosleep(&req, (struct timespec *)NULL);
}

void wait(ThreadHandle aThreadHandle)
{
	pthread_join(aThreadHandle->thread, 0);
}

void release(ThreadHandle aThreadHandle)
{
	delete aThreadHandle;
}

int getTimeMillis()
{
	struct timespec spec;
	clock_gettime(CLOCK_REALTIME, &spec);
	return spec.tv_sec * 1000 + (int)(spec.tv_nsec / 1.0e6);
}
#endif

#endif // __cplusplus >= 201703L

static void poolWorker(void *aParam)
{
	Pool *myPool = (Pool *)aParam;
	while (myPool->mRunning)
	{
		PoolTask *t = myPool->getWork();
		if (!t)
		{
			sleep(1);
		}
		else
		{
			t->work();
		}
	}
}

Pool::Pool()
{
	mRunning = 0;
	mThreadCount = 0;
	mThread = nullptr;
	mWorkMutex = nullptr;
	mRobin = 0;
	mMaxTask = 0;
	for (int i = 0; i < MAX_THREADPOOL_TASKS; i++)
		mTaskArray[i] = nullptr;
}

Pool::~Pool()
{
	mRunning = 0;
	int i;
	for (i = 0; i < mThreadCount; i++)
	{
		wait(mThread[i]);
		release(mThread[i]);
	}
	delete[] mThread;
	if (mWorkMutex)
		destroyMutex(mWorkMutex);
}

void Pool::init(int aThreadCount)
{
	if (aThreadCount > 0)
	{
		mMaxTask = 0;
		mWorkMutex = createMutex();
		mRunning = 1;
		mThreadCount = aThreadCount;
		mThread = new ThreadHandle[aThreadCount];
		int i;
		for (i = 0; i < mThreadCount; i++)
		{
			mThread[i] = createThread(poolWorker, this);
		}
	}
}

void Pool::addWork(PoolTask *aTask)
{
	if (mThreadCount == 0)
	{
		aTask->work();
	}
	else
	{
		if (mWorkMutex)
			lockMutex(mWorkMutex);
		if (mMaxTask == MAX_THREADPOOL_TASKS)
		{
			// If we're at max tasks, do the task on calling thread
			// (we're in trouble anyway, might as well slow down adding more work)
			if (mWorkMutex)
				unlockMutex(mWorkMutex);
			aTask->work();
		}
		else
		{
			mTaskArray[mMaxTask] = aTask;
			mMaxTask++;
			if (mWorkMutex)
				unlockMutex(mWorkMutex);
		}
	}
}

PoolTask *Pool::getWork()
{
	PoolTask *t = nullptr;
	if (mWorkMutex)
		lockMutex(mWorkMutex);
	if (mMaxTask > 0)
	{
		int r = mRobin % mMaxTask;
		mRobin++;
		t = mTaskArray[r];
		mTaskArray[r] = mTaskArray[mMaxTask - 1];
		mMaxTask--;
	}
	if (mWorkMutex)
		unlockMutex(mWorkMutex);
	return t;
}
} // namespace Thread
} // namespace SoLoud
