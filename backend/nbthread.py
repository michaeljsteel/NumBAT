# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import threading
import multiprocessing
import traceback
import queue
import time
import datetime
#import numbat


class CalcThread(threading.Thread):
    """Runs the calculation function f_work(task) in a separate thread on one or more tasks
    popped from queue.Queue q_work with results pushed to queue q_result.

    First element of the task tuple drawn from the queue should be
    an integer(8) identifying which piece of work is being done.

    An additional optional queue q_work_noshare contains any tasks that may involve
    non-thread-safe tasks, such as matplotlib plotting. Only one thread should be supplied
    a non-empty q_work_noshare. Obviously, this queue should be be small compared to the
    main work queue or there will be minimal parallel advantage.

    """

    def __init__(self, q_work, q_result, f_work, q_work_noshare=None, verbose=False):
        """Optional q_work_noshare for non-thread-safe tasks.
        Other threads get an empty such queue and will only draw tasks from q_work."""
        self.q_work = q_work
        if q_work_noshare is None:
            self.q_work_noshare = (
                queue.Queue()
            )  # Non Thread-1 workers get their own empty queue
        else:
            self.q_work_noshare = q_work_noshare
        self.q_result = q_result

        self.f_work = f_work

        self.td = 0
        self.ptd = 0

        self.verbose = verbose

        # only one thread may call matplotlib
        # self.can_plot = can_plot
        self.doing_plot_work = False

        threading.Thread.__init__(self)

    def run(self):
        while True:
            task = None
            try:
                if self.verbose:
                    print(
                        "{0} encountered remaining work queue lengths: ({1}, {2}).".format(
                            self.name, self.q_work_noshare.qsize(), self.q_work.qsize()
                        )
                    )

                if (
                    not self.q_work_noshare.empty()
                ):  # I should exhaust q_work_noshare first
                    task = (
                        self.q_work_noshare.get()
                    )  # not empty and sole task_func so must succeed
                    self.doing_plot_work = True
                else:
                    task = self.q_work.get_nowait()

            except queue.Empty:
                if self.verbose:
                    print(f"{self.name} out of work, wrapping up.")
                break

            if self.verbose:
                print(
                    f"{self.name} is starting task i={task[0]}, needs plotting: ",
                    self.doing_plot_work,
                )

            try:
                MAGIC_TASK = (
                    -10,
                    -10,
                    -10,
                    -10,
                )  # trigger an easy first task to get around GIL issues
                if (
                    task == MAGIC_TASK
                ):  # complete immediately to try and release the GIL
                    print(f"{self.name} is doing the magic task!")
                    # res = MAGIC_TASK

                else:
                    res = self.f_work(task)  # Here is the main piece of work

                    if self.verbose:
                        print(f"{self.name} produced outcome: ", res)

                    self.q_result.put(res)

                if self.doing_plot_work:
                    self.q_work_noshare.task_done()
                    self.doing_plot_work = False
                    self.ptd += 1
                    # print('\n\n{0} has called ptd {1} times.'.format(self.name, self.ptd))
                else:
                    self.q_work.task_done()
                    self.td += 1
                    # print('\n\n{0} has called td {1} times.'.format(self.name, self.td))

            except Exception as err:
                print(
                    f"\n\n{self.name} has encountered an exception on task {task[0]}:",
                    str(err),
                )
                traceback.print_exception(err)
                break

        return


def run_in_current_thread(task_func, q_result, q_work, verbose=False):
    '''Simply executes task_func(tsk) on the tasks in q_work in a single
    non-threaded loop, placing the results in q_result.'''

    total_tasks = q_work.qsize()
    for i_task in range(total_tasks):
        print('\n\n----------------------------------------')
        print(f'Starting task {i_task+1} of {total_tasks}.\n\n')
        task = q_work.get(block=True, timeout=5)
        res = task_func(task)  # Here is the main piece of work
        q_result.put(res)


def launch_worker_threads_and_wait(num_workers, task_func, q_result,
                                   q_work, q_work_noshare=None, verbose=False):

    if num_workers < 1:  # avoid separate thread if num_workers <= 0
        run_in_current_thread(task_func, q_result, q_work, verbose)
        return

    report_progress = True
    # verbose=False

    total_tasks = q_work.qsize()
    if q_work_noshare is not None:
        total_tasks += q_work_noshare.qsize()
    # Launch threads and keep copies

    num_workers = min(num_workers, total_tasks)
    if report_progress:
        print(f"Assigning {total_tasks} tasks across {num_workers} threads.")

    threads = []
    for i in range(num_workers):
        print("Looking at starting thread", i)

        # If there _is_ a non-Null q_work_noshare, only thread 1 should see it
        qwns = q_work_noshare if (i == 0) else None

        th = CalcThread(q_work, q_result, task_func, qwns, verbose)
        print("Starting thread:", th.name)
        th.start()
        print("Done starting thread:", th.name)
        threads.append(th)

    tm_st = time.time()

    # Wait for all work to be done

    if report_progress:
        pause = 5  # progress reporting interval in seconds
        pause_max = 120  # maximum reporting interval in seconds
    else:
        pause = None

    for th in threads:  # wait on worker threads.
        while True:
            th.join(pause)
            if report_progress:  # Wake up now and again to report progress
                if (
                    pause < pause_max
                ):  # Early on, we report frequently, then less often over time
                    pause *= 1.25

                tasks_done = total_tasks - q_work.qsize()
                if q_work_noshare is not None:
                    tasks_done -= q_work_noshare.qsize()
                tm_cu = time.time()
                frac_done = tasks_done / total_tasks
                if frac_done > 0:
                    tm_togo = (tm_cu - tm_st) * (1 - frac_done) / frac_done
                    tm_togo_s = datetime.timedelta(seconds=round(tm_togo))

                    print(
                        f"\nTasks completed: {tasks_done}/{total_tasks} = {frac_done*100:.1f} %.",
                        f"  Estimated time remaining: {tm_togo_s}.",
                    )
                else:
                    print(f"(\nTasks completed: 0/{total_tasks}.")

            if not th.is_alive():
                print(
                    f"Main thread has joined thread {th.name}. ",
                    f"There are {threading.active_count()-1} workers remaining.",
                )
                break

    print("All tasks completed and worker threads terminated.")


class CalcProcess(multiprocessing.Process):
    """Runs the calculation function f_work(task) in a separate process on one or more tasks
    popped from multiprocessing.JoinableQueue q_work with results pushed to queue q_result.

    First element of the task tuple drawn from the queue should be
    an integer(8) identifying which piece of work is being done.

    An additional optional queue q_work_noshare contains any tasks that may involve
    non-thread-safe tasks, such as matplotlib plotting. Only one thread should be supplied
    a non-empty q_work_noshare. Obviously, this queue should be be small compared to the
    main work queue or there will be minimal parallel advantage.

    """

    def __init__(self, q_work, q_result, f_work, verbose=False):

        self.q_work = q_work
        self.q_result = q_result

        self.f_work = f_work

        self.td = 0
        self.ptd = 0

        self.verbose = verbose

        multiprocessing.Process.__init__(self)

    def run(self):
        while True:
            task = None
            try:
                if self.verbose:
                    print(
                        f"{self.name} encountered remaining work queue length: {self.q_work.qsize()}."
                    )

                # This timed wait to get a chance at the queue is inelegant
                # Did have a get_nowait(), but then when all processes start at once,
                #  some miss out on getting access and decide to give up.
                # Better would be that nothing quits until it actually measures a qsize()=0
                #  or else, we wait to acquire a lock on the queue, before testing get_nowait()
                task = self.q_work.get(block=True, timeout=5)

            except queue.Empty:
                if self.verbose:
                    print(f"{self.name} out of work, wrapping up.")
                break

            if self.verbose:
                print(f"{self.name} is starting task i={task[0]}.")

            try:
                res = self.f_work(task)  # Here is the main piece of work

                # if self.verbose: print('{0} produced outcome:'.format(self.name), res)

                self.q_result.put(res)

                self.q_work.task_done()
                self.td += 1

            except Exception as err:
                print(
                    f"\n\n{self.name} has encountered an exception on task {task[0]}:",
                    str(err),
                )
                traceback.print_exception(err)
                break

        print(f"{self.name} is running off the end")


def launch_worker_processes_and_wait(num_workers, task_func,
                                     q_result, q_work, verbose=False):

    if num_workers < 1:  # avoid separate thread if num_workers <= 0
        run_in_current_thread(task_func, q_result, q_work, verbose)
        return

    report_progress = True
    # verbose=True

    # TODO: turned off access to numbat to avoid circulat import dependency
    #do_multiproc = numbat.NumBATApp().can_multiprocess()


    # Launch processes and keep copies

    if do_multiproc:
        total_tasks = q_work.qsize()
        num_workers = min(num_workers, total_tasks)

        if report_progress:
            print(f"Assigning {total_tasks} tasks across {num_workers} processes.")
    else:
        print("Performing calculation in single processor mode.")
        num_workers = 1

    processes = []
    for _ in range(num_workers):
        pr = CalcProcess(q_work, q_result, task_func, verbose)
        print("Starting process:", pr.name)
        if do_multiproc:
            pr.start()
            processes.append(pr)
        else:
            pr.run()

    tm_st = time.time()

    # Wait for all work to be done

    if report_progress:
        pause = 5  # progress reporting interval in seconds
        pause_max = 120  # maximum reporting interval in seconds
    else:
        pause = None

    for pr in processes:  # wait on worker threads.
        while True:
            pr.join(pause)
            if report_progress:  # Wake up now and again to report progress
                if (
                    pause < pause_max
                ):  # Early on, we report frequently, then less often over time
                    pause *= 1.25

                tm_cu = time.time()

                tasks_started = total_tasks - q_work.qsize()
                tasks_completed = q_result.qsize()
                frac_started = tasks_started / total_tasks
                frac_completed = tasks_completed / total_tasks

                print(
                    f"\nTasks commenced: {tasks_started}/{total_tasks} = {frac_started*100:.1f} %;",
                    f" completed: {tasks_completed}/{total_tasks} = {frac_completed*100:.1f} %;",
                    end="",
                )

                if frac_completed > 0:
                    tm_togo = (tm_cu - tm_st) * (1 - frac_completed) / frac_completed
                    tm_togo_s = datetime.timedelta(seconds=round(tm_togo))
                    print(f"  Estimated time remaining: {tm_togo_s}.")
                else:
                    print("\n")

            if not pr.is_alive():
                print(
                    f"Main process has joined process {pr.name}. ",
                    f"There are {len(multiprocessing.active_children())} workers remaining.",
                )
                break

    print("All tasks completed and worker process terminated.")
