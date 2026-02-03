from dataclasses import dataclass
from typing import Literal, Union

@dataclass
class FlowRankConfig:
    """
        Configuration for FlowRank algorithm.
        q- int: Number of nearest neighbors to consider for constructing the k-NN graph.
        r- int: Number of neighbors to consider during the ascending random walk process.

    """


    q: int = 20
    r: int = 10


@dataclass
class StableCoreConfig:

    """
        densification:
        - "k-nn":  local k-nearest-neighbor densification
        - "rw":    random-walk-based densification
        - False:   skip refinement

    """

    auto_select_core_frac: bool = False
    core_frac: float = 0.15
    resolution: float = 1.0
    ng_num: int = 10
    densification: Union[Literal["k-nn", "rw"], bool] = False


@dataclass
class FineGrainedConfig:

    stable_core_densification: Union[Literal["k-nn", "rw"], bool] = False
    auto_select_core_frac: bool = False
    core_frac: float=0.2
    fine_grain_densification : Union[Literal["k-nn", "rw"], bool] = False
    starting_resolution: float = 1.5


@dataclass
class ClusterConfig:

    """
    Configuration for clustering core points using the Leiden algorithm.
    resolution: float
        Resolution parameter for Leiden clustering.
    """
    densification: Union[Literal["k-nn", "rw"], bool] = False
    resolution: float = 1.0
    auto_select_resolution: bool = False


@dataclass
class PropagationConfig:

    """
    Configuration for label propagation from core points.
    normalized: bool
        Whether to use a normalized graph Laplacian

    connectivity: {"in", "out", "sym"}
        - "in": Use incoming edges for propagation.
        - "out": Use outgoing edges for propagation.
        - "sym": Use symmetrized edges for propagation.


    mode: {"default", "momentum"}
        - "default": Standard label propagation.
        - "momentum": Momentum-based propagation for potentially faster convergence.

    alpha: float
        Damping factor for propagation, typically between 0 and 1.

    tol: float
        Tolerance for convergence.

    max_iter: int
        Maximum number of iterations for the propagation to perform.
    """


    normalized: bool = True
    connectivity: Literal["in", "out", "sym"] = "sym"
    mode : Literal["default", "momentum"] = "momentum"
    alpha: float = 0.1 #Slowly progressing walk to account for imbalanced cluster sizes.
    tol: float = 1e-6
    max_iter: int = 50


from dataclasses import dataclass, field, replace
from typing import Union, Literal, Dict


from dataclasses import dataclass, field, replace
from typing import Union, Literal, Dict
from .config import (
    FlowRankConfig,
    StableCoreConfig,
    FineGrainedConfig,
    ClusterConfig,
    PropagationConfig,
)


@dataclass
class CoreSpectConfig:
    """
    Unified configuration for the CoreSpect pipeline.

    Global parameters
    -----------------
    densify : {"k-nn", "rw", False}, default=False
        If not False, enables densification of all relevant stages
        ("k-nn" or "rw" options).
    q, r : int
        Passed to FlowRankConfig.
    core_frac : float
        Passed to both StableCoreConfig and FineGrainedConfig.
    auto_select_core_frac : bool
        Passed to both StableCoreConfig and FineGrainedConfig.
    granularity : float
        Sets StableCoreConfig.resolution and FineGrainedConfig.starting_resolution.
    resolution : float
        Sets ClusterConfig.resolution.
    auto_select_resolution : bool
        Sets ClusterConfig.auto_select_resolution.
    """

    # Stage-specific configs (internal)
    flowrank: FlowRankConfig = field(default_factory=FlowRankConfig)
    stable: StableCoreConfig = field(default_factory=StableCoreConfig)
    fine: FineGrainedConfig = field(default_factory=FineGrainedConfig)
    cluster: ClusterConfig = field(default_factory=ClusterConfig)
    propagation: PropagationConfig = field(default_factory=PropagationConfig)

    # --- Unified interface (user-facing global flags) ---
    densify: Union[Literal["k-nn", "rw"], bool] = False

    q: int = None
    r: int = None
    core_frac: float = None
    auto_select_core_frac: bool = None
    granularity: float = None
    resolution: float = None
    auto_select_resolution: bool = None

    # ------------------------------------------------------------------
    def configure(self):
        """Apply global flags across all relevant stage configs."""
        # --- 1. Densification ---
        if self.densify in ("k-nn", "rw"):
            self.stable = replace(self.stable, densification=self.densify)
            self.fine = replace(
                self.fine,
                stable_core_densification=self.densify,
                fine_grain_densification=self.densify,
            )
            self.cluster = replace(self.cluster, densification=self.densify)
        elif self.densify is False:
            self.stable = replace(self.stable, densification=False)
            self.fine = replace(
                self.fine,
                stable_core_densification=False,
                fine_grain_densification=False,
            )
            self.cluster = replace(self.cluster, densification=False)
        else:
            raise ValueError("densify must be one of {'k-nn', 'rw', False}")

        # --- 2. FlowRank parameters ---
        updates = {}
        if self.q is not None:
            updates["q"] = self.q
        if self.r is not None:
            updates["r"] = self.r
        if updates:
            self.flowrank = replace(self.flowrank, **updates)

        # --- 3. Core fraction controls ---
        updates = {}
        if self.core_frac is not None:
            updates["core_frac"] = self.core_frac
        if self.auto_select_core_frac is not None:
            updates["auto_select_core_frac"] = self.auto_select_core_frac
        if updates:
            self.stable = replace(self.stable, **updates)
            self.fine = replace(self.fine, **updates)

        # --- 4. Granularity mapping ---
        if self.granularity is not None:
            self.stable = replace(self.stable, resolution=self.granularity)
            self.fine = replace(self.fine, starting_resolution=self.granularity)

        # --- 5. Cluster resolution controls ---
        updates = {}
        if self.resolution is not None:
            updates["resolution"] = self.resolution
        if self.auto_select_resolution is not None:
            updates["auto_select_resolution"] = self.auto_select_resolution
        if updates:
            self.cluster = replace(self.cluster, **updates)

        return self

    # ------------------------------------------------------------------
    def unpack(self) -> Dict[str, object]:
        """Return configs in a format directly usable by CorespectModel."""
        return dict(
            flowrank_cfg=self.flowrank,
            stable_cfg=self.stable,
            fine_core_cfg=self.fine,
            cluster_cfg=self.cluster,
            prop_cfg=self.propagation,
        )

    # ------------------------------------------------------------------
    def summary(self):
        """Human-readable summary."""
        print("CoreSpect configuration")
        print(f"  Global densification: {self.densify}")
        print(f"  q={self.q}, r={self.r}, core_frac={self.core_frac}, granularity={self.granularity}")
        print(f"  resolution={self.resolution}, auto_select_resolution={self.auto_select_resolution}")
        print("  Stage configs:")
        for name, sub in [
            ("FlowRank", self.flowrank),
            ("StableCore", self.stable),
            ("FineGrained", self.fine),
            ("Cluster", self.cluster),
            ("Propagation", self.propagation),
        ]:
            print(f"   - {name}: {vars(sub)}")
