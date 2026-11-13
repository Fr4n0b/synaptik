"""
Adaptive Dose Escalation Engines: Traditional 3+3, BOIN, and mTPI.

Implements rule-based 3+3 and Bayesian model-assisted designs (Bayesian Optimal
Interval [BOIN] and modified Toxicity Probability Interval [mTPI]) for First-in-Human
intrathecal clinical trials in rare neurogenetic disorders.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, List, Optional


class EscalationAction(str, Enum):
    """Recommended clinical trial action for subsequent patient cohort."""
    ESCALATE = "escalate"
    RETAIN = "retain"
    DE_ESCALATE = "de_escalate"
    STOP_TRIAL = "stop_trial"
    DECLARE_MTD = "declare_mtd"


class EscalationDesign(str, Enum):
    """Dose escalation algorithmic framework."""
    THREE_PLUS_THREE = "3+3"
    BOIN = "boin"
    MTPI = "mtpi"


@dataclass
class CohortOutcome:
    """Historical or observed outcome of an enrolled cohort."""
    dose_level: int
    dose_mg: float
    n_enrolled: int
    n_dlts: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "dose_level": self.dose_level,
            "dose_mg": self.dose_mg,
            "n_enrolled": self.n_enrolled,
            "n_dlts": self.n_dlts,
            "observed_dlt_rate": round(self.n_dlts / self.n_enrolled, 4) if self.n_enrolled > 0 else 0.0,
        }


@dataclass
class EscalationDecision:
    """Dose escalation recommendation and Bayesian safety assessment."""
    current_dose_level: int
    current_dose_mg: float
    action: EscalationAction
    next_dose_level: Optional[int]
    next_dose_mg: Optional[float]
    mtd_declared_dose_level: Optional[int]
    mtd_declared_dose_mg: Optional[float]
    observed_toxicity_rate: float
    posterior_overdose_prob: float
    is_dose_eliminated: bool
    eliminated_dose_levels: List[int]
    rationale: str
    design: str = "boin"
    target_phi: float = 0.25
    lambda_e: float = 0.197
    lambda_d: float = 0.298
    trial_stopped: bool = False
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    @property
    def posterior_prob(self) -> float:
        """Alias for posterior_overdose_prob."""
        return self.posterior_overdose_prob

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "design": self.design,
            "target_phi": self.target_phi,
            "lambda_e": round(self.lambda_e, 4),
            "lambda_d": round(self.lambda_d, 4),
            "current_dose_level": self.current_dose_level,
            "current_dose_mg": self.current_dose_mg,
            "recommendation": self.action.value,
            "recommended_dose_level": self.next_dose_level,
            "recommended_dose_mg": self.next_dose_mg,
            "mtd_declared_dose_level": self.mtd_declared_dose_level,
            "mtd_declared_dose_mg": self.mtd_declared_dose_mg,
            "observed_dlt_rate": round(self.observed_toxicity_rate, 4),
            "posterior_excess_toxicity_prob": round(self.posterior_overdose_prob, 4),
            "is_dose_eliminated": self.is_dose_eliminated,
            "eliminated_dose_levels": self.eliminated_dose_levels,
            "is_trial_stopped": self.trial_stopped,
            "safety_stopping_triggered": self.trial_stopped or self.is_dose_eliminated,
            "rationale": self.rationale,
        }


# ---------------------------------------------------------------------------
# Pure-Python Special Functions (Incomplete Beta Function for Posterior Prob)
# ---------------------------------------------------------------------------


def _betacf(a: float, b: float, x: float, max_iter: int = 200, eps: float = 1e-12) -> float:
    """Continued fraction evaluation of the regularized incomplete Beta function."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        # Even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c

        # Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        del_val = d * c
        h *= del_val
        if abs(del_val - 1.0) < eps:
            break
    return h


def betai(a: float, b: float, x: float) -> float:
    """
    Pure-Python regularized incomplete Beta function I_x(a, b).
    Calculates integral_0^x t^(a-1) (1-t)^(b-1) dt / B(a, b) with 0 external dependencies.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    else:
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def calculate_beta_overdose_prob(dlts: int, n: int, target_phi: float, alpha_prior: float = 1.0, beta_prior: float = 1.0) -> float:
    """
    Calculate Bayesian posterior probability of overdosing: P(p > phi | y, n).
    Under conjugate Beta(alpha_prior, beta_prior), posterior is Beta(alpha_prior + y, beta_prior + n - y).
    """
    a = alpha_prior + dlts
    b = beta_prior + n - dlts
    # P(p > phi) = 1.0 - I_phi(a, b)
    p_cum = betai(a, b, target_phi)
    prob = 1.0 - p_cum
    return max(0.0, min(1.0, prob))


# ---------------------------------------------------------------------------
# BOIN Escalation Engine (Liu & Yuan 2015)
# ---------------------------------------------------------------------------


class BOINEscalationEngine:
    """
    Bayesian Optimal Interval (BOIN) dose-escalation engine.

    Calculates optimal likelihood escalation and de-escalation boundaries
    (lambda_e, lambda_d) based on target toxicity rate phi, and enforces
    Bayesian conjugate posterior dose elimination rules P(p > phi | data) > 0.95.
    """

    def __init__(
        self,
        target_phi: float = 0.25,
        dose_levels: Optional[List[float]] = None,
        safety_cutoff: float = 0.95,
    ) -> None:
        if target_phi <= 0.0 or target_phi >= 1.0:
            raise ValueError("Target toxicity rate phi must be strictly between 0 and 1")
        if safety_cutoff <= 0.0 or safety_cutoff >= 1.0:
            raise ValueError("Safety cutoff must be strictly between 0 and 1")

        self.target_phi = float(target_phi)
        self.safety_cutoff = float(safety_cutoff)
        self.dose_levels = list(dose_levels) if dose_levels else [10.0, 25.0, 50.0, 100.0, 150.0]

        # Standard BOIN boundary derivations:
        # phi_1 = 0.6 * phi, phi_2 = 1.4 * phi
        phi_1 = 0.6 * self.target_phi
        phi_2 = 1.4 * self.target_phi

        self.lambda_e = math.log((1.0 - phi_1) / (1.0 - self.target_phi)) / math.log(
            (self.target_phi * (1.0 - phi_1)) / (phi_1 * (1.0 - self.target_phi))
        )
        self.lambda_d = math.log((1.0 - self.target_phi) / (1.0 - phi_2)) / math.log(
            (phi_2 * (1.0 - self.target_phi)) / (self.target_phi * (1.0 - phi_2))
        )

    def evaluate_cohort(
        self,
        dose_level: int,
        n: int,
        dlts: int,
        history: Optional[List[CohortOutcome]] = None,
    ) -> EscalationDecision:
        """
        Evaluate cohort observation at dose_level and determine next action.

        Parameters
        ----------
        dose_level : int
            1-indexed dose level tier (1 .. len(dose_levels)).
        n : int
            Total evaluable patients at this dose level.
        dlts : int
            Number of dose-limiting toxicities observed at this level.
        history : list of CohortOutcome, optional
            Cumulative multi-cohort trial history.
        """
        if dose_level < 1 or dose_level > len(self.dose_levels):
            raise ValueError(f"dose_level must be between 1 and {len(self.dose_levels)}")
        if n <= 0:
            raise ValueError("Cohort size n must be positive")
        if dlts < 0 or dlts > n:
            raise ValueError("dlts must be between 0 and n")

        p_hat = dlts / n
        current_dose_mg = self.dose_levels[dose_level - 1]
        p_over = calculate_beta_overdose_prob(dlts, n, self.target_phi)

        eliminated_levels: List[int] = []
        is_eliminated = False
        trial_stopped = False

        # Safety stopping rule: if P(p_j > phi) > 0.95 with n >= 3
        if n >= 3 and p_over > self.safety_cutoff:
            is_eliminated = True
            eliminated_levels = list(range(dose_level, len(self.dose_levels) + 1))
            if dose_level == 1:
                trial_stopped = True
                return EscalationDecision(
                    current_dose_level=1,
                    current_dose_mg=current_dose_mg,
                    action=EscalationAction.STOP_TRIAL,
                    next_dose_level=None,
                    next_dose_mg=None,
                    mtd_declared_dose_level=None,
                    mtd_declared_dose_mg=None,
                    observed_toxicity_rate=p_hat,
                    posterior_overdose_prob=p_over,
                    is_dose_eliminated=True,
                    eliminated_dose_levels=eliminated_levels,
                    rationale=(
                        f"Excessive toxicity observed at starting dose Level 1: {dlts}/{n} DLTs "
                        f"(Posterior P(p > {self.target_phi:.2f}) = {p_over:.4f} > {self.safety_cutoff:.2f}). "
                        f"Starting dose level eliminated; trial terminated immediately."
                    ),
                    design="boin",
                    target_phi=self.target_phi,
                    lambda_e=self.lambda_e,
                    lambda_d=self.lambda_d,
                    trial_stopped=True,
                )
            else:
                next_lvl = dose_level - 1
                return EscalationDecision(
                    current_dose_level=dose_level,
                    current_dose_mg=current_dose_mg,
                    action=EscalationAction.DE_ESCALATE,
                    next_dose_level=next_lvl,
                    next_dose_mg=self.dose_levels[next_lvl - 1],
                    mtd_declared_dose_level=None,
                    mtd_declared_dose_mg=None,
                    observed_toxicity_rate=p_hat,
                    posterior_overdose_prob=p_over,
                    is_dose_eliminated=True,
                    eliminated_dose_levels=eliminated_levels,
                    rationale=(
                        f"Dose Level {dose_level} eliminated due to excessive toxicity "
                        f"(Posterior P(p > {self.target_phi:.2f}) = {p_over:.4f} > {self.safety_cutoff:.2f}). "
                        f"De-escalating to Level {next_lvl} ({self.dose_levels[next_lvl - 1]} mg)."
                    ),
                    design="boin",
                    target_phi=self.target_phi,
                    lambda_e=self.lambda_e,
                    lambda_d=self.lambda_d,
                    trial_stopped=False,
                )

        # Standard BOIN interval decision rules:
        if p_hat <= self.lambda_e:
            if dose_level == len(self.dose_levels):
                action = EscalationAction.RETAIN
                next_lvl = dose_level
                rationale = (
                    f"Observed DLT rate {dlts}/{n} ({p_hat:.3f}) <= lambda_e ({self.lambda_e:.3f}). "
                    f"Already at maximum planned dose Level {dose_level} ({current_dose_mg} mg); retain current level."
                )
            else:
                action = EscalationAction.ESCALATE
                next_lvl = dose_level + 1
                rationale = (
                    f"Observed DLT rate {dlts}/{n} ({p_hat:.3f}) <= lambda_e ({self.lambda_e:.3f}). "
                    f"Escalate to dose Level {next_lvl} ({self.dose_levels[next_lvl - 1]} mg)."
                )
        elif self.lambda_e < p_hat <= self.lambda_d:
            action = EscalationAction.RETAIN
            next_lvl = dose_level
            rationale = (
                f"Observed DLT rate {dlts}/{n} ({p_hat:.3f}) is within interval "
                f"({self.lambda_e:.3f}, {self.lambda_d:.3f}]. Retain current dose Level {dose_level} ({current_dose_mg} mg)."
            )
        else:
            # p_hat > self.lambda_d
            if dose_level == 1:
                action = EscalationAction.STOP_TRIAL
                next_lvl = None
                trial_stopped = True
                rationale = (
                    f"Observed DLT rate {dlts}/{n} ({p_hat:.3f}) > lambda_d ({self.lambda_d:.3f}) at lowest dose Level 1. "
                    f"De-escalation impossible; trial stopped."
                )
            else:
                action = EscalationAction.DE_ESCALATE
                next_lvl = dose_level - 1
                rationale = (
                    f"Observed DLT rate {dlts}/{n} ({p_hat:.3f}) > lambda_d ({self.lambda_d:.3f}). "
                    f"De-escalate to dose Level {next_lvl} ({self.dose_levels[next_lvl - 1]} mg)."
                )

        next_mg = self.dose_levels[next_lvl - 1] if next_lvl is not None else None

        return EscalationDecision(
            current_dose_level=dose_level,
            current_dose_mg=current_dose_mg,
            action=action,
            next_dose_level=next_lvl,
            next_dose_mg=next_mg,
            mtd_declared_dose_level=None,
            mtd_declared_dose_mg=None,
            observed_toxicity_rate=p_hat,
            posterior_overdose_prob=p_over,
            is_dose_eliminated=is_eliminated,
            eliminated_dose_levels=eliminated_levels,
            rationale=rationale,
            design="boin",
            target_phi=self.target_phi,
            lambda_e=self.lambda_e,
            lambda_d=self.lambda_d,
            trial_stopped=trial_stopped,
        )


# ---------------------------------------------------------------------------
# Traditional 3+3 Escalation Engine
# ---------------------------------------------------------------------------


class ThreePlusThreeEngine:
    """
    Deterministic rule-based traditional 3+3 dose escalation engine.

    Standard protocol:
    - 0/3 DLTs: Escalate to next dose level.
    - 1/3 DLTs: Expand current cohort by 3 additional patients (N=6).
    - 1/6 DLTs: Escalate to next dose level.
    - >=2/3 or >=2/6 DLTs: MTD exceeded. Halt escalation; declare prior dose MTD.
    """

    def __init__(self, dose_levels: Optional[List[float]] = None) -> None:
        self.dose_levels = list(dose_levels) if dose_levels else [10.0, 25.0, 50.0, 100.0, 150.0]

    def evaluate_cohort(
        self,
        dose_level: int,
        n: int,
        dlts: int,
        history: Optional[List[CohortOutcome]] = None,
    ) -> EscalationDecision:
        """Evaluate 3+3 decision tree based on observed DLTs and cohort size."""
        if dose_level < 1 or dose_level > len(self.dose_levels):
            raise ValueError(f"dose_level must be between 1 and {len(self.dose_levels)}")
        if n not in (3, 6):
            if n <= 0:
                raise ValueError("Cohort size n must be positive")
            # We allow arbitrary n for robustness, but standard 3+3 evaluates at 3 and 6
        if dlts < 0 or dlts > n:
            raise ValueError("dlts must be between 0 and n")

        p_hat = dlts / n
        current_dose_mg = self.dose_levels[dose_level - 1]
        p_over = calculate_beta_overdose_prob(dlts, n, 0.25)

        action: EscalationAction
        next_lvl: Optional[int]
        mtd_lvl: Optional[int] = None
        mtd_mg: Optional[float] = None
        trial_stopped = False
        rationale: str

        if n == 3:
            if dlts == 0:
                if dose_level == len(self.dose_levels):
                    action = EscalationAction.DECLARE_MTD
                    next_lvl = None
                    mtd_lvl = dose_level
                    mtd_mg = current_dose_mg
                    trial_stopped = True
                    rationale = f"0/3 DLTs observed at maximum planned dose Level {dose_level}. Declare Level {dose_level} as MTD; trial completed."
                else:
                    action = EscalationAction.ESCALATE
                    next_lvl = dose_level + 1
                    rationale = f"0/3 DLTs observed at dose Level {dose_level}. Escalate to Level {next_lvl} ({self.dose_levels[next_lvl - 1]} mg)."
            elif dlts == 1:
                action = EscalationAction.RETAIN
                next_lvl = dose_level
                rationale = f"1/3 DLTs observed at dose Level {dose_level}. Expand cohort to N=6 at current dose Level {dose_level} ({current_dose_mg} mg)."
            else:
                # dlts >= 2
                action = EscalationAction.STOP_TRIAL
                next_lvl = None
                trial_stopped = True
                if dose_level == 1:
                    mtd_lvl = None
                    mtd_mg = None
                    rationale = f">=2/3 DLTs ({dlts}/3) observed at starting dose Level 1. Tolerance exceeded; MTD not identified. Trial terminated."
                else:
                    mtd_lvl = dose_level - 1
                    mtd_mg = self.dose_levels[mtd_lvl - 1]
                    rationale = f">=2/3 DLTs ({dlts}/3) observed at dose Level {dose_level}. MTD exceeded. Declare prior Level {mtd_lvl} ({mtd_mg} mg) as MTD."
        else:
            # n == 6 or expanded cohort
            if dlts <= 1:
                if dose_level == len(self.dose_levels):
                    action = EscalationAction.DECLARE_MTD
                    next_lvl = None
                    mtd_lvl = dose_level
                    mtd_mg = current_dose_mg
                    trial_stopped = True
                    rationale = f"{dlts}/6 DLTs observed at maximum dose Level {dose_level}. Acceptable tolerance; declare Level {dose_level} as MTD."
                else:
                    action = EscalationAction.ESCALATE
                    next_lvl = dose_level + 1
                    rationale = f"{dlts}/6 DLTs observed at dose Level {dose_level}. Acceptable tolerance; escalate to Level {next_lvl} ({self.dose_levels[next_lvl - 1]} mg)."
            else:
                # dlts >= 2
                action = EscalationAction.STOP_TRIAL
                next_lvl = None
                trial_stopped = True
                if dose_level == 1:
                    mtd_lvl = None
                    mtd_mg = None
                    rationale = f">=2/6 DLTs ({dlts}/{n}) observed at starting dose Level 1. Tolerance exceeded; MTD not identified. Trial terminated."
                else:
                    mtd_lvl = dose_level - 1
                    mtd_mg = self.dose_levels[mtd_lvl - 1]
                    rationale = f">=2/6 DLTs ({dlts}/{n}) observed at dose Level {dose_level}. MTD exceeded. Declare prior Level {mtd_lvl} ({mtd_mg} mg) as MTD."

        next_mg = self.dose_levels[next_lvl - 1] if next_lvl is not None else None

        return EscalationDecision(
            current_dose_level=dose_level,
            current_dose_mg=current_dose_mg,
            action=action,
            next_dose_level=next_lvl,
            next_dose_mg=next_mg,
            mtd_declared_dose_level=mtd_lvl,
            mtd_declared_dose_mg=mtd_mg,
            observed_toxicity_rate=p_hat,
            posterior_overdose_prob=p_over,
            is_dose_eliminated=(action == EscalationAction.STOP_TRIAL),
            eliminated_dose_levels=list(range(dose_level, len(self.dose_levels) + 1)) if action == EscalationAction.STOP_TRIAL else [],
            rationale=rationale,
            design="3+3",
            target_phi=0.25,
            lambda_e=0.17,
            lambda_d=0.33,
            trial_stopped=trial_stopped,
        )


# ---------------------------------------------------------------------------
# Modified Toxicity Probability Interval (mTPI) Engine
# ---------------------------------------------------------------------------


class MTPIEngine:
    """
    Modified Toxicity Probability Interval (mTPI) dose escalation engine.

    Partitions the toxicity space [0, 1] into three intervals:
    - Underdosing: I_E = [0, phi - eps)
    - Target:      I_T = [phi - eps, phi + eps]
    - Overdosing:  I_D = (phi + eps, 1]
    and selects the action with the maximal Unit Probability Mass (UPM).
    """

    def __init__(
        self,
        target_phi: float = 0.25,
        epsilon: float = 0.05,
        dose_levels: Optional[List[float]] = None,
        safety_cutoff: float = 0.95,
    ) -> None:
        if target_phi <= 0.0 or target_phi >= 1.0:
            raise ValueError("Target toxicity rate phi must be strictly between 0 and 1")
        if epsilon <= 0.0 or (target_phi - epsilon) <= 0.0 or (target_phi + epsilon) >= 1.0:
            raise ValueError("Epsilon produces degenerate toxicity intervals")

        self.target_phi = float(target_phi)
        self.epsilon = float(epsilon)
        self.safety_cutoff = float(safety_cutoff)
        self.dose_levels = list(dose_levels) if dose_levels else [10.0, 25.0, 50.0, 100.0, 150.0]

    def evaluate_cohort(
        self,
        dose_level: int,
        n: int,
        dlts: int,
        history: Optional[List[CohortOutcome]] = None,
    ) -> EscalationDecision:
        """Evaluate mTPI decision based on Unit Probability Mass (UPM)."""
        if dose_level < 1 or dose_level > len(self.dose_levels):
            raise ValueError(f"dose_level must be between 1 and {len(self.dose_levels)}")
        if n <= 0:
            raise ValueError("Cohort size n must be positive")
        if dlts < 0 or dlts > n:
            raise ValueError("dlts must be between 0 and n")

        p_hat = dlts / n
        current_dose_mg = self.dose_levels[dose_level - 1]
        p_over = calculate_beta_overdose_prob(dlts, n, self.target_phi)

        # Safety stopping rule
        if n >= 3 and p_over > self.safety_cutoff:
            eliminated = list(range(dose_level, len(self.dose_levels) + 1))
            trial_stopped = (dose_level == 1)
            next_lvl = None if trial_stopped else (dose_level - 1)
            return EscalationDecision(
                current_dose_level=dose_level,
                current_dose_mg=current_dose_mg,
                action=EscalationAction.STOP_TRIAL if trial_stopped else EscalationAction.DE_ESCALATE,
                next_dose_level=next_lvl,
                next_dose_mg=self.dose_levels[next_lvl - 1] if next_lvl else None,
                mtd_declared_dose_level=None,
                mtd_declared_dose_mg=None,
                observed_toxicity_rate=p_hat,
                posterior_overdose_prob=p_over,
                is_dose_eliminated=True,
                eliminated_dose_levels=eliminated,
                rationale=f"mTPI safety stopping: P(p > {self.target_phi}) = {p_over:.4f} > {self.safety_cutoff}.",
                design="mtpi",
                target_phi=self.target_phi,
                trial_stopped=trial_stopped,
            )

        # Partition intervals
        low = self.target_phi - self.epsilon
        high = self.target_phi + self.epsilon
        a = 1.0 + dlts
        b = 1.0 + n - dlts

        prob_e = betai(a, b, low)
        prob_t = betai(a, b, high) - betai(a, b, low)
        prob_d = 1.0 - betai(a, b, high)

        upm_e = prob_e / low
        upm_t = prob_t / (high - low)
        upm_d = prob_d / (1.0 - high)

        if upm_e >= upm_t and upm_e >= upm_d:
            if dose_level == len(self.dose_levels):
                action = EscalationAction.RETAIN
                next_lvl = dose_level
                rationale = f"mTPI UPM_E ({upm_e:.3f}) is highest at maximum Level {dose_level}; retain current dose."
            else:
                action = EscalationAction.ESCALATE
                next_lvl = dose_level + 1
                rationale = f"mTPI UPM_E ({upm_e:.3f}) is highest (underdosing); escalate to Level {next_lvl}."
        elif upm_t >= upm_e and upm_t >= upm_d:
            action = EscalationAction.RETAIN
            next_lvl = dose_level
            rationale = f"mTPI UPM_T ({upm_t:.3f}) is highest (target toxicity reached); retain Level {dose_level}."
        else:
            if dose_level == 1:
                action = EscalationAction.STOP_TRIAL
                next_lvl = None
                rationale = f"mTPI UPM_D ({upm_d:.3f}) is highest at lowest Level 1; stop trial."
            else:
                action = EscalationAction.DE_ESCALATE
                next_lvl = dose_level - 1
                rationale = f"mTPI UPM_D ({upm_d:.3f}) is highest (overdosing); de-escalate to Level {next_lvl}."

        next_mg = self.dose_levels[next_lvl - 1] if next_lvl else None

        return EscalationDecision(
            current_dose_level=dose_level,
            current_dose_mg=current_dose_mg,
            action=action,
            next_dose_level=next_lvl,
            next_dose_mg=next_mg,
            mtd_declared_dose_level=None,
            mtd_declared_dose_mg=None,
            observed_toxicity_rate=p_hat,
            posterior_overdose_prob=p_over,
            is_dose_eliminated=False,
            eliminated_dose_levels=[],
            rationale=rationale,
            design="mtpi",
            target_phi=self.target_phi,
            trial_stopped=(action == EscalationAction.STOP_TRIAL),
        )
