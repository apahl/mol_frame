import time
import pickle

try:
    import pandas as pd

    PANDAS = True
except ImportError:
    PANDAS = False


def is_interactive_ipython():
    try:
        get_ipython()  # type: ignore
        ipy = True
        print("> interactive IPython session.")
    except NameError:
        ipy = False
    return ipy


IPYTHON = is_interactive_ipython()


if IPYTHON:
    from IPython.core.display import HTML, Javascript, display
    import uuid


def format_seconds(seconds):
    seconds = int(seconds)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    t_str = "{:02.0f}h {:02d}m {:02d}s".format(h, m, s)
    return t_str


class ProgCtr:
    """A ProgressCounter Class"""

    def __init__(self, val=0):
        self.val = val

    def __call__(self):
        return self.val

    def inc(self, x=1):
        self.val += 1


class Result:
    def __init__(self, headers=["Result", "Value"]):
        self.headers = headers
        self.list = []

    def __str__(self):
        col0_max = max([len(x[0]) for x in self.list])
        col1_max = max([len(x[1]) for x in self.list])
        if col0_max < len(self.headers[0]):
            col0_max = len(self.headers[0])
        if col1_max < len(self.headers[1]):
            col1_max = len(self.headers[1])
        sep = "―" * (col0_max + col1_max + 2)
        out = [f"{self.headers[0]:{col0_max}s}  {self.headers[1]:>{col1_max}s}"]
        out.append(sep)
        for r in self.list:
            out.append(f"{r[0]:{col0_max}s}  {r[1]:>{col1_max}s}")
        return "\n".join(out)

    def __repr__(self):
        return self.__str__()

    def to_html(self):
        if PANDAS and IPYTHON:
            return pd.DataFrame.from_records(self.list, columns=self.headers)
        else:
            return ""

    def add(self, *res, show=True):
        """Add one or more result tuples to the instance.
        Show the instance."""
        for r in res:
            if isinstance(r, str):
                if len(self.list) > 0:
                    self.list.append((" ", " "))
                self.list.append((r, " "))
                continue
            r = list(r)
            if len(r) < 2:
                r.append(" ")
                if len(self.list) > 0:
                    self.list.append((" ", " "))
            else:
                r[0] = "• " + r[0]
                if isinstance(r[1], int):
                    r[1] = str(r[1])
                elif isinstance(r[1], float):
                    r[1] = f"{r[1]:.3f}"
            self.list.append(r)
        if show:
            print(self.__str__())

    def remove(self, n):
        """Remove the n last entries."""
        self.list = self.list[:-n]

    def tmp(self, *res):
        """Temporarily add an entry.
        During development. Add, show, delete."""
        self.add(res)
        self.remove(len(res))

    def clear(self):
        self.list = []


class Progressbar:
    """A class to display a Javascript progressbar in the IPython notebook."""

    def __init__(self, ctr=0, end=100, color="#43ace8"):
        self.ctr = ctr
        self.end = end
        self.bar_id = str(uuid.uuid4())
        self.eta_id = str(uuid.uuid4())
        # possible colours: #94CAEF (blue from HTML reports),
        #                   #d6d2d0 (grey from window decorations)
        #
        self.pb = HTML(
            """
            <table style="border: none;"><tbody><tr style="border: none;">
            <td style="border: none;"><div style="border: 1px solid black; height:6px; width:500px">
                <div id="{}" style="background-color:{}; height:4px; width:0%">&nbsp;</div>
            </div></td><td style="border: none;">&nbsp;ETA:&nbsp;</td><td style="border: none;"><div id="{}" width=100px></div></td>
            </tr></tbody></table>
            """.format(
                self.bar_id, color, self.eta_id
            )
        )
        self.prev_time = 0.0
        self.start_time = time.time()
        display(self.pb)

    def update(self, perc, force=False):
        """update the progressbar
        in: progress in percent"""
        # make sure that the update function is not called too often:
        self.cur_time = time.time()
        if force or (self.cur_time - self.prev_time >= 1.0):
            if perc > 100:
                perc = 100
            if perc >= 25:
                eta = (100 - perc) * (self.cur_time - self.start_time) / perc
                eta_str = format_seconds(eta)
            else:
                eta_str = "..."
            self.prev_time = self.cur_time
            display(
                Javascript(
                    """
                                    $('div#{}').width('{}%');
                                    $('div#{}').text('{}');
                            """.format(
                        self.bar_id, perc, self.eta_id, eta_str
                    )
                )
            )

    def inc(self):
        self.ctr += 1
        self.update(100 * self.ctr / self.end)

    def done(self):
        """finalize with a full progressbar for aesthetics"""
        display(
            Javascript(
                """
                                $('div#{}').width('{}%');
                                $('div#{}').text('{}');
                            """.format(
                    self.bar_id, 100, self.eta_id, "done"
                )
            )
        )


def listify(s, sep=None, as_int=True):
    """A helper func for the Jupyter Notebook,
    which generates a correctly formatted list out of pasted text."""
    to_number = int if as_int else float
    result = []
    if sep is None:
        if "\n" in s[:-2]:
            sep = "\n"
        elif "," in s[:-2]:
            sep = ","
        elif "\t" in s[:-2]:
            sep = "\t"
        else:
            sep = " "
    if s.startswith("["):
        s = s[1:]
    if s.endswith("]"):
        s = s[:-1]
    lst = s.split(sep)
    for el in lst:
        if len(el) == 0:
            continue
        try:
            el = to_number(el)
        except ValueError:
            pass
        result.append(el)
    return result


def print_list(lst):
    for el in lst:
        print(el)


def save_object(obj, fn="object.pkl"):
    """Save a generic python object through pickling."""
    with open(fn, "wb") as f:
        pickle.dump(obj, f)


def load_object(fn="object.pkl"):
    """Load and return a previously pickled object."""
    with open(fn, "rb") as f:
        obj = pickle.load(f)
    return obj
