from enum import Enum


class Trace(Enum):
    UP = 1
    LEFT = 2
    DIAG = 3
    NONE = 4
    

class TracedScore:
    def __init__(self, t: Trace, s: int):
        self.trace = t
        self.score = s

    def __repr__(self) -> str:
        if self.trace == Trace.UP:
            a = r"\uparrow"
        elif self.trace == Trace.LEFT:
            a = r"\leftarrow"
        elif self.trace == Trace.DIAG:
            a = r"\nwarrow"
        else:
            a = ""

        return f"({a} {self.score})"
