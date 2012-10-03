"""
================
Warping programs
================

The ``@command`` decorator makes binding command line executables
easy to write and easy to use.

To wrap a program, write a function that takes whatever arguments
you will need to run the program like if it where on the shell.
The function should return a dictionary containing several keys:

    * mandatory: ``arguments``
    * optional: ``stdin``
    * optional: ``return_value``.

Firstly, ``arguments`` should be a list of strings which is
the actual command and arguments to be executed (e.g. ``["touch", filename]``).

Secondly, ``stdin`` should be a value to feed the subprocess once it is launched.

Thirdly, ``return_value`` be a value to return, or a callable
object which takes a stdout plus strerr as parameters and returns the value
that will be passed back to the user when this program is run.
You can also simply specify ``"stdout"`` to have the output of the
process returned directly.

For example, to wrap ``touch``, we write a one argument function that
takes the filename of the file to touch, and apply the ``@command``
decorator to it::

    from plumbing import command

    @command
    def touch(filename):
        return {"arguments": ["touch", filename],
                "return_value": filename}

We can now call this function directly::

     f = touch("myfile")

The value returned by touch is ``"myfile"``, the name of
the touched file.

A more complicated example would include binding the BLASTP algorithm::

    from plumbing import command

    @command
    def blastp(database, sequences, **kwargs):
        \"\"\"Will launch the 'blastp' algorithm using the NCBI executable.

       :param database: The path to the database to blast against.
       :param sequences: A fasta formated string.
       :param **kwargs: Extra parameters that will be passed to the executable
                        For instance, you could specify "e=1e-20".
       :returns: A list of top hits in blast format.
       \"\"\"
        return {"arguments": ["blastall", "-p", "blastp", "-d" database] +
                             [a for k,v in kwargs.items() for a in ('-'+k,v)],
                "stdin": sequences)
                "return_value": 'stdout'}

Often you want to call a function, but not block when it returns
so you can run several in parallel. ``@command`` also creates a
method ``parallel`` which does this. The return value is a
Future object with a single method: ``wait()``. When you call
``wait()``, it blocks until the program finishes, then returns the
same value that you would get from calling the function directly.
So, to touch two files, and not block until both commands have
started, you would write::

    a = blastp.parallel("nr", open("fasta1").read(), e=1e-20)
    b = blastp.parallel("nr", open("fasta2").read(), e=1e-20)
    hitsA = a.wait()
    hitsB = b.wait()

The ``parallel`` method will runs processes without blocking.
"""

# Built-in modules #
import sys, subprocess, time

# Variables #
PARRALEL_JOBS = []

################################################################################
def start_process(args):
    """Run a process using subprocess module"""
    try: return subprocess.Popen(args, bufsize=-1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    except OSError: raise Exception("Program '%s' does not seem to exist in your $PATH." % args[0])

################################################################################
def pause_for_parralel_jobs(update_interval=2):
    """Wait until all parallel jobs are done and print a status update"""
    global PARRALEL_JOBS
    while True:
        PARRALEL_JOBS = [job for job in PARRALEL_JOBS if not job.finished]
        if not PARRALEL_JOBS:
            sys.stdout.write("\r\033[K")
            sys.stdout.flush()
            return
        sys.stdout.write("\r    %i parallel jobs still running.\033[K" % len(PARRALEL_JOBS))
        sys.stdout.flush()
        time.sleep(update_interval)

################################################################################
class CommandFailed(Exception):
    """Thrown when a program bound by ``@command``
    exits with a value other than ``0``."""

    def __init__(self, args, stderr=None):
        message = "Running '%s' failed." % " ".join(args)
        if stderr: message += " The error reported is:\n\n" + stderr
        Exception.__init__(self, message)

################################################################################
class command(object):
    """Decorator used to wrap external programs."""

    def __init__(self, function):
        self.function = function
        self.__doc__ = function.__doc__
        self.__name__ = function.__name__

    def __call__(self, *args, **kwargs):
        """Run a program locally, and block until it completes."""
        # Call the user defined function #
        cmd_dict = self.function(*args, **kwargs)
        args = cmd_dict['arguments']
        # Start a process #
        proc = start_process(args)
        # Wait until completion #
        try: stdout, stderr = proc.communicate(cmd_dict.get("stdin"))
        except KeyboardInterrupt as err:
            print "You aborted the process pid %i. It was: %s " % (proc.pid, args)
            raise err
        # Check for failure #
        if proc.returncode != 0: raise CommandFailed(args, stderr)
        # Return result #
        result = cmd_dict.get("return_value")
        if callable(result): return result(stdout, stderr)
        elif result == 'stdout': return stdout
        else: return result

    def parallel(self, *args, **kwargs):
        """Run a program and return a Future object."""
        # Call the user defined function #
        cmd_dict = self.function(*args, **kwargs)
        # Start a process #
        proc = start_process(cmd_dict['arguments'])
        # Write the standard in #
        if 'stdin' in cmd_dict:
            proc.stdin.write(cmd_dict["stdin"])
            proc.stdin.close()
        # The Future object takes it from here #
        future = Future(proc, cmd_dict)
        # Let's keep a reference of it #
        PARRALEL_JOBS.append(future)
        # Hand it back to the user #
        return future

################################################################################
class Future(object):
    """Object returned when functions decorated with ``@command``
    are executed in parallel with ``parallel()``."""

    def __init__(self, proc, cmd_dict):
        self.proc = proc
        self.cmd_dict = cmd_dict

    @property
    def finished(self):
        if self.proc.poll() == None: return False
        else: return True

    def wait(self):
        # Wait until completion #
        try: return_code = self.proc.wait()
        except KeyboardInterrupt as err:
            print "You aborted the process pid %i. It was: %s " % (self.proc.pid, self.cmd_dict["arguments"])
            raise err
        # Read result #
        stdout, stderr = self.proc.stdout.read(),  self.proc.stderr.read()
        # Check for failure #
        if return_code != 0: raise CommandFailed(self.cmd_dict["arguments"], stderr)
        # Return result #
        result = self.cmd_dict.get("return_value")
        if callable(result): return result(stdout, stderr)
        elif result == 'stdout': return stdout
        else: return result