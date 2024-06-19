"""Implementation of the Needleman-Wunsch algorithm
with incremental latex matrix output support.
"""

from enum import Flag, auto
from typing import Optional


class Trace(Flag):
    UP = 0b0001
    LEFT = 0b0010
    DIAG = 0b0100
    NONE = 0b1000


class TracedScore:
    def __init__(self, t: Trace, s: int):
        self.trace = t
        self.score = s

    def __lt__(self, other):
        return self.score < other.score

    def __le__(self, other):
        return self.score <= other.score

    def __gt__(self, other):
        return self.score > other.score

    def __ge__(self, other):
        return self.score >= other.score

    def __eq__(self, other):
        return self.score == other.score

    def __ne__(self, other):
        return self.score != other.score

    def __repr__(self) -> str:
        return f"TracedScore({self.trace}, {self.score})"

    def __str__(self) -> str:
        """Formats self for use as LaTeX matrix element"""
        diag = r"\nwarrow" if self.trace & Trace.DIAG else r"\phantom\nwarrow"
        up = r"\uparrow" if self.trace & Trace.UP else r"\phantom\uparrow"
        left = r"\leftarrow" if self.trace & Trace.LEFT else r"\phantom\leftarrow"
        neg = r"\llap{^-}" if self.score < 0 else ""

        return f"^{{{diag} {up}}}_{{{left} {neg}{abs(self.score)}}}"


class NeedlemanWunsch:
    def __init__(
        self, a: str, b: str, match: int = 1, mismatch: int = -1, gap_penalty: int = 2
    ):
        if gap_penalty <= 0:
            raise ValueError(f"gap_penalty must be positive. Got {gap_penalty}.")

        self.a = f"–{a.replace('-', '')}"
        self.b = f"–{b.replace('-', '')}"
        self.match = match
        self.mismatch = mismatch
        self.gap_penalty = abs(gap_penalty)

        self._row: int = 0
        self._col: int = 0
        self._matrix: list[list[TracedScore]] = [[] for _ in range(len(self.a))]

        self._tracedback: bool = False
        self.traceback_indices: list[tuple[int, int]] = []

        self.i = None
        self.j = None

    def traceback(self) -> tuple[str, str]:
        """Performs the traceback, storing the traceback indices in self.traceback_indices

        Returns
        -------
        a_aln: str
            Aligned version of self.a
        b_aln: str
            Aligned version of self.b
        """
        if self.i < len(self.a) - 1 or self.j < len(self.b) - 1:
            raise ValueError("Need to iterate before calling self.traceback()")

        while self.i > 0 and self.j > 0:
            self.traceback_indices.append((self.i, self.j))
            if self.trace() & Trace.DIAG:
                self.i -= 1
                self.j -= 1
            elif self.trace() & Trace.LEFT:
                self.j -= 1
            elif self.trace() & Trace.UP:
                self.i -= 1
            else:
                raise RuntimeError(f"i={self.i}, j={self.j}")

        self._tracedback = True

        a_aln = []
        b_aln = []
        i_prev, j_prev = 0, 0
        for i, j in reversed(self.traceback_indices):
            if i == i_prev + 1 and j == j_prev + 1:
                a_aln.append(self.a[i])
                b_aln.append(self.b[j])
            elif i == i_prev + 1 and j == j_prev:
                a_aln.append(self.a[i])
                b_aln.append("-")
            elif i == i_prev and j == j_prev + 1:
                a_aln.append("-")
                b_aln.append(self.b[j])
            else:
                RuntimeError("This line should be unreachable.")

            i_prev, j_prev = i, j

        return "".join(a_aln), "".join(b_aln)

    def score(self, i: Optional[int] = None, j: Optional[int] = None) -> int:
        """Get the score at position [i, j] in the matrix

        Parameters
        ----------
        i: int
            Row index (by default use current index)
        j: int
            Column index (by default use current index)

        Returns
        -------
        score: int

        Raises
        ------
        IndexError
            Care must be taken to ensure the relevant value has already been calulated
        """
        if i is None:
            i = self.i

        if j is None:
            j = self.j

        return self._matrix[i][j].score

    def trace(self, i: Optional[int] = None, j: Optional[int] = None) -> Trace:
        """Get the trace at position [i, j] in the matrix

        Parameters
        ----------
        i: int
            Row index (by default use current index)
        j: int
            Column index (by default use current index)

        Returns
        -------
        trace: Trace

        Raises
        ------
        IndexError
            Care must be taken to ensure the relevant value has already been calulated
        """
        if i is None:
            i = self.i

        if j is None:
            j = self.j

        return self._matrix[i][j].trace

    def __repr__(self) -> str:
        return (
            f"NeedlemanWunsch({self.a}, {self.b}, match={self.match}, "
            f"mismatch={self.mismatch}, gap_penalty={self.gap_penalty})"
        )

    def __str__(self) -> str:
        """Formats current state of algorithm as a LaTeX matrix"""
        colours = (("black", "orange"), ("red", "blue"))
        if not self._tracedback:
            rows = (
                f"\\text{{{c}}} & {r' & '.join(f"\\color{{{colours[self.i - i][self.j - j]}}}{{{s}}}" if (self.i > 0 and self.j > 0) and (0 <= self.i - i <= 1) and (0 <= self.j - j <= 1) else str(s) for j, s in enumerate(r))}"
                for i, (c, r) in enumerate(zip(self.a, self._matrix))
            )
        else:
            rows = (
                f"\\text{{{c}}} & {r' & '.join(f"\\boxed{{{s}}}" if (i, j) in self.traceback_indices else str(s) for j, s in enumerate(r))}"
                for i, (c, r) in enumerate(zip(self.a, self._matrix))
            )

        matrix = "\\\\\n".join(rows)
        return f"$$\\begin{{matrix}}\n & \\text{{{'} & \\text{'.join(self.b)}}}\\\\\n{matrix}\n\\end{{matrix}}$$"

    def initialise_scores(self) -> bool:
        """Initialises first row and column scores and sets pointer.
        Does nothing if already initialised.

        Returns
        -------
        initialised: bool
            True if this is the first time self.initialise_scores has been called
        """
        if self._row:
            return False

        self._matrix = [[TracedScore(Trace.NONE, 0)]]
        self._matrix[0].extend(
            TracedScore(Trace.LEFT, s)
            for s in range(self.mismatch, len(self.b) * self.mismatch, self.mismatch)
        )
        self._matrix.extend(
            [TracedScore(Trace.UP, s)]
            for s in range(self.mismatch, len(self.a) * self.mismatch, self.mismatch)
        )
        self.i = 1
        self.j = 0
        self._row = 1
        self._col = 1

        return True

    def __iter__(self):
        return self

    def __next__(self):
        if self._row >= len(self.a):
            raise StopIteration

        if self.initialise_scores():
            return TracedScore(Trace.NONE, 0)

        match_mismatch = (
            self.match if self.a[self._row] == self.b[self._col] else self.mismatch
        )

        traced_scores = sorted(
            (
                TracedScore(
                    Trace.LEFT, self._matrix[self._row][-1].score - self.gap_penalty
                ),
                TracedScore(
                    Trace.UP,
                    self._matrix[self._row - 1][self._col].score - self.gap_penalty,
                ),
                TracedScore(
                    Trace.DIAG,
                    self._matrix[self._row - 1][self._col - 1].score + match_mismatch,
                ),
            )
        )

        best_score = traced_scores.pop()
        for traced_score in traced_scores:
            if traced_score.score == best_score.score:
                best_score.trace |= traced_score.trace

        self._matrix[self._row].append(best_score)

        self.i = self._row
        self.j = self._col

        self._col += 1
        if self._col >= len(self.b):
            self._col = 1
            self._row += 1

        return best_score
