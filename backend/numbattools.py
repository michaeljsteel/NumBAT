import threading
import multiprocessing
import traceback
import queue
import time
import datetime
import scipy.integrate as sciint

import numpy as np


import numbat
import reporting


def almost_zero(x, tol=1e-10):
    return abs(x) < tol


def almost_unity(x, tol=1e-10):
    return abs(1 - x) < tol


def np_min_max(v):
    return np.min(v), np.max(v)


def int2d(mat):
    return sciint.simpson(sciint.simpson(mat))


def process_fortran_return(resm, msg):

    fort_err, fort_mesg = resm[-2:]
    if fort_err:
        fort_mesg = str(fort_mesg, "utf-8")  # fort_mesg comes back as a byte string.
        reporting.report_and_exit(
            f"Fortran error in {msg}: \n"
            f" NumBAT Fortran error code = {fort_err}. \n Message: \n {fort_mesg}"
        )
        return None  # dummy
    else:  # everything is fine
        return resm[:-2]


class CalcThread(threading.Thread):
    """Runs the calculation function f_work(task) in a separate thread on one or more tasks
    popped from queue.Queue q_work with results pushed to queue q_result.

    First element of the task tuple drawn from the queue should be
    an integer identifying which piece of work is being done.

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
                    )  # not empty and sole caller so must succeed
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


def launch_worker_threads_and_wait(
    num_threads, caller, q_result, q_work, q_work_noshare=None, verbose=False
):

    # TODO: avoid separate thread if num_threads = 1
    report_progress = True
    # verbose=False

    total_tasks = q_work.qsize()
    if q_work_noshare is not None:
        total_tasks += q_work_noshare.qsize()
    # Launch threads and keep copies

    num_threads = min(num_threads, total_tasks)
    if report_progress:
        print(f"Assigning {total_tasks} tasks across {num_threads} threads.")

    threads = []
    for i in range(num_threads):
        print("Looking at starting thread", i)

        # If there _is_ a non-Null q_work_noshare, only thread 1 should see it
        qwns = q_work_noshare if (i == 0) else None

        th = CalcThread(q_work, q_result, caller, qwns, verbose)
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
    an integer identifying which piece of work is being done.

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


def launch_worker_processes_and_wait(
    num_processes, caller, q_result, q_work, verbose=False
):

    # TODO: avoid separate thread if num_processes = 1
    report_progress = True
    # verbose=True

    do_multiproc = numbat.NumBATApp().can_multiprocess()

    # Launch processes and keep copies

    if do_multiproc:
        total_tasks = q_work.qsize()
        num_processes = min(num_processes, total_tasks)

        if report_progress:
            print(f"Assigning {total_tasks} tasks across {num_processes} processes.")
    else:
        print("Performing calculation in single processor mode.")
        num_processes = 1

    processes = []
    for _ in range(num_processes):
        pr = CalcProcess(q_work, q_result, caller, verbose)
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
