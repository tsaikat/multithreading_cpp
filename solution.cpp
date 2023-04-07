#ifndef __PROGTEST__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <array>
#include <iterator>
#include <set>
#include <list>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <deque>
#include <memory>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <stdexcept>
#include <condition_variable>
#include <pthread.h>
#include <semaphore.h>
#include "progtest_solver.h"
#include "sample_tester.h"
using namespace std;
#endif /* __PROGTEST__ */
// -------------------------------------------------------------------------------------------------------------------------------------------------------------

#define BUFFER_SIZE 50

class ProblemPackWrapper
{
public:
    ProblemPackWrapper(AProblemPack problems, ACompany company, int order)
        : m_problems(problems), m_company(company), m_order(order)
    {
        unsolved = m_problems->m_Problems.size();
    }

    AProblemPack m_problems;
    ACompany m_company;
    int m_order;
    int unsolved;
    std::mutex mtx;
};

using ProblemPackPlus = std::shared_ptr<ProblemPackWrapper>;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

class ProblemPackPlusComparator
{
public:
    bool operator()(const ProblemPackPlus &a, const ProblemPackPlus &b) const

    {
        return a->m_order > b->m_order;
    }
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
class ProblemWrapper
{
public:
    ProblemWrapper(AProblem theproblem, ProblemPackPlus pack, bool lastProblem = false)
        : problem(theproblem), origin(pack), last(lastProblem)
    {
    }

    AProblem problem;
    ProblemPackPlus origin;
    bool last;
};

using ProblemPlus = std::shared_ptr<ProblemWrapper>;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

class CProblemQueue
{
public:
    void push(ProblemPlus problems)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        cv_full.wait(lock, [this]()
                     { return (m_queue.size() < BUFFER_SIZE); });
        m_queue.push(problems);

        //std::cout << "+ Pushed" << std::endl;

        cv_empty.notify_one();
    }

    ProblemPlus pop()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        cv_empty.wait(lock, [this]()
                      { //std::cout << "Sleeping..." << std::endl;
                    return !m_queue.empty() ; });

        ProblemPlus problems = m_queue.front();
        m_queue.pop();

        //std::cout << "- Popped" << std::endl;
        
        cv_full.notify_one();
        return problems;
    }

    bool empty()
    {
        return m_queue.empty();
    }

    std::queue<ProblemPlus> m_queue;
    std::mutex m_mutex;
    std::condition_variable cv_full, cv_empty;
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

class COptimizer
{
public:
    static bool usingProgtestSolver(void)
    {
        return true;
    }

    static void checkAlgorithm(AProblem problem)
    {
        // dummy implementation if usingProgtestSolver() returns true
    }

    void start(int threadCount)
    {
        runningWorker = threadCount;
        runningProd = m_companies.size();

        for (int i = 0; i < threadCount; ++i)
        {
            m_threads.emplace_back(&COptimizer::worker, this);
        }
        intiateCommunication();
    }

    void stop(void)
    {
        //std::cout << "STOP STARTED" << std::endl;
        for (auto &thread : m_threads)
        {
            if (thread.joinable())
                thread.join();
        }
        //std::cout << "STOPPED and SUCCESS!" << std::endl;
    }

    void addCompany(ACompany company)
    {
        m_companies.push_back(company);
    }

    void intiateCommunication()
    {
        for (auto company : m_companies)
        {
            m_threads.emplace_back(&COptimizer::receiveProblemPacks, this, company);
            m_threads.emplace_back(&COptimizer::sendSolvedPacks, this, company);
        }
    }

private:
    void receiveProblemPacks(ACompany company)
    {
        int order = 0;
        while (true)
        {
            AProblemPack basicPorblemPack = company->waitForPack();

            if (basicPorblemPack)
            {
                ProblemPackPlus problemPack = std::make_shared<ProblemPackWrapper>(basicPorblemPack, company, order++);
                auto& rawProblem = problemPack->m_problems->m_Problems;
                for ( auto& p : rawProblem ) {
                    ProblemPlus wrappedProblem =std::make_shared<ProblemWrapper> (p, problemPack);  
                    m_problemQueue.push(wrappedProblem);
                }
            }
            else
            {
                // thread dies when nothing else to receive from company
                unique_lock<mutex> lock(runningMtx);
                int activeProd = runningProd;
                int activeWorkers = runningWorker;
                lock.unlock();

                if (activeProd == 1)
                {
                    for (int i = 0; i < activeWorkers; i++)
                    {
                        // creating poison so the worker threads can die
                        ProblemPlus wrappedProblem =std::make_shared<ProblemWrapper> (nullptr, nullptr, true);  
                        m_problemQueue.push(wrappedProblem);
                    }
                }
                unique_lock<mutex> lock2(runningMtx);
                runningProd--;
                lock2.unlock();
                
                break;
            }
        }
        std::cout << "Receive problems: ENDED" << std::endl;
    }




    void sendSolvedPacks(ACompany company)
    {
        std::priority_queue<ProblemPackPlus, std::vector<ProblemPackPlus>,
                            ProblemPackPlusComparator>
            pQueue;

        int next_order = 0;
        // unique_lock<mutex> ul(runningMtx);
        while (true)
        {   
            // ul.unlock();
            std::unique_lock<std::mutex> lock(m_solvedPacksMutex);
            m_cv.wait(lock, [this, &company]
                      { 
                        //std::cout << "Sleeping...insolvedPack" << runningWorker << std::endl;
                        return m_solvedPacks[company].size() > 0 || (runningWorker == 0) ;});
            //std::cout << "up...insolvedPack" << std::endl;
            if (runningWorker == 0 && m_solvedPacks[company].empty()) {
                std::cout << "-->" << std::endl;
                break;
            }

            std::vector<ProblemPackPlus> buffer = std::move(m_solvedPacks[company]);
            lock.unlock();
            m_cv.notify_one();

            for (auto &pack : buffer)
            {
                pQueue.push(pack);
            }

            while (pQueue.top()->m_order == next_order)
            {
                company->solvedPack(pQueue.top()->m_problems);
                std::cout << "SOLVED---> " <<pQueue.top()->m_order << std::endl;
                pQueue.pop();
                ++next_order;
            }
            // ul.lock();
        }
        
        std::cout << "Send SOLVED: ENDED" << std::endl;
    }

    void deliverToSolveThread(ProblemPlus problem)
    {
        std::unique_lock<mutex> originLock(problem->origin->mtx);
        int unsolved = --(problem->origin->unsolved);
        originLock.unlock();
        
        if ( unsolved == 0 ) { // all solved in this pack
            std::unique_lock<std::mutex> lock(m_solvedPacksMutex);
            m_solvedPacks[problem->origin->m_company].push_back(problem->origin);
            std::cout << "Delivered to solve"<< "------> " << problem->origin->m_order << std::endl;
            m_cv.notify_all();
        }
    }





    void worker()
    {
        AProgtestSolver solver = createProgtestSolver();
        
        if (!solver->hasFreeCapacity()) {
            unique_lock<mutex> workerVar(runningMtx);
            runningWorker--;
            if (runningWorker == 0) m_cv.notify_all();
            workerVar.unlock();
            return;
        }

        std::vector<ProblemPlus> toBeSolved;

        unique_lock<mutex> lock(runningMtx);
        while (!(m_problemQueue.empty() && runningProd == 0))
        {
            lock.unlock();
            ProblemPlus problem = m_problemQueue.pop();

            if (problem->last)
            {
                solver->solve();
                for (auto& item : toBeSolved) {
                    deliverToSolveThread(item);
                }
                // toBeSolved.clear(); no need to clear cause we're exiting
                break;
            }

            solver->addProblem(problem->problem);
            toBeSolved.push_back(problem);

            if (!solver->hasFreeCapacity()) {
                //don't have capacity, handle it
                solver->solve();
                for (auto& item : toBeSolved) {
                    deliverToSolveThread(item);
                }
                toBeSolved.clear();
                solver = createProgtestSolver();
                if (!solver->hasFreeCapacity()) break;
            }
            
            lock.lock();
        }

        unique_lock<mutex> workerVar(runningMtx);
        runningWorker--;
        if (runningWorker == 0) m_cv.notify_all();
        workerVar.unlock();
        std::cout << "Consumer : end" << std::endl;
    }

    std::vector<ACompany> m_companies;
    CProblemQueue m_problemQueue;
    std::vector<std::thread> m_threads;

    int runningProd;
    int runningWorker;
    std::mutex runningMtx;


    std::mutex m_solvedPacksMutex;
    std::condition_variable m_cv;
    std::map<ACompany, std::vector<ProblemPackPlus>> m_solvedPacks;
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef __PROGTEST__
int main(void)
{
    COptimizer optimizer;
    ACompanyTest company = std::make_shared<CCompanyTest>();
    ACompanyTest company2 = std::make_shared<CCompanyTest>();
    optimizer.addCompany(company);
    optimizer.addCompany(company2);
    optimizer.start(34);
    optimizer.stop();
    if (!company->allProcessed())
        throw std::logic_error("(some) problems were not correctly processsed");
    if (company->allProcessed())
        std::cout << "all problems were processsed" << std::endl;
    return 0;
}
#endif /* __PROGTEST__ */
