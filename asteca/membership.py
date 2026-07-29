from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .cluster import Cluster

import numpy as np


class Membership:
    """Define a :py:class:`Membership` object.

    This object is used as a container for membership probabilities methods. Currently
    two methods are included:

    :py:meth:`bayesian` : The algorithm was described in detail in
    the `article <https://doi.org/10.1051/0004-6361/201424946>`__
    were we originally introduced ``ASteCA``. The method requires ``(RA, DEC)``
    data and will use any extra data dimensions stored in the
    :py:class:`Cluster <asteca.cluster.Cluster>` object, i.e.: photometry,
    proper motions, and parallax.

    :py:meth:`fastmp` : The algorithm was described in detail in
    the `article <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__
    were we introduced the `Unified Cluster Catalogue (UCC) <https://ucc.ar/>`__.
    The method requires proper motions, and parallax data dimensions stored in the
    :py:class:`Cluster <asteca.cluster.Cluster>` object. Photometry data is not
    employed.

    :param my_field: :py:class:`Cluster <asteca.cluster.Cluster>` object with the
        loaded data for the observed field
    :type my_field: Cluster
    """

    def __init__(
        self,
        my_field: Cluster,
    ) -> None:
        self.my_field = my_field

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.my_field.verbose > level:
            print(mssg)

    def bayesian(
        self,
        N_runs: int = 1000,
        eq_to_gal: bool = False,
        seed: int | None = None,
    ) -> None:
        """Assign membership probabilities.

        Estimate the probability of being a true cluster member for all observed
        stars, using a Bayesian algorithm. The ``radec_c`` and ``radius`` attributes
        are required to be present in the
        :py:class:`Cluster <asteca.cluster.Cluster>` object.

        :param N_runs: Maximum number of runs
        :type N_runs: int
        :param eq_to_gal: Convert ``(RA, DEC)`` to ``(lon, lat)``. Useful for clusters
            with large ``DEC`` values to reduce the frame's distortion
        :type eq_to_gal: bool
        :param seed: Random seed. If ``None`` a random integer will be generated and used
        :type seed: int | None

        :raises AttributeError: If either required attributes
            are missing from the :py:class:`Cluster <asteca.cluster.Cluster>` object
        :raises ValueError: If the ``N_runs`` parameter is less than 10
        """
        self._vp("\nRunning Bayesian DA")

        # These two are checked separately to avoid type checking issues
        if self.my_field.ra is None or self.my_field.dec is None:
            raise AttributeError("Missing required 'cluster' attributes: (ra, dec)")
        required_attrs = ("N_cluster", "radec_c", "radius")
        missing = [attr for attr in required_attrs if not hasattr(self.my_field, attr)]
        if missing:
            raise AttributeError(
                f"Missing required 'cluster' attributes: {', '.join(missing)}"
            )

        if N_runs < 10:
            raise ValueError("Parameter 'N_runs' should be > 10")

        # Set seed
        self.rng = np.random.default_rng(seed)

        from .modules import cluster_priv as cp
        from .modules.bayesian_da import bayesian_mp

        xv, yv = self.my_field.ra, self.my_field.dec
        xy_center = self.my_field.radec_c
        cent_str = "radec_c "
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            # Add center as last element
            xv_t = np.array(list(xv) + [xy_center[0]])
            yv_t = np.array(list(yv) + [xy_center[1]])
            # Convert
            xv_t, yv_t = cp.radec2lonlat(xv_t, yv_t)
            # Extract center value
            xv, yv = xv_t[:-1], yv_t[:-1]
            xy_center = (xv_t[-1], yv_t[-1])
            cent_str = "lonlat_c"

        # Generate input data array
        X = [xv, yv]
        e_X = []
        if self.my_field.mag is not None:
            X.append(self.my_field.mag)
            e_X.append(self.my_field.e_mag)
        if self.my_field.color is not None:
            X.append(self.my_field.color)
            e_X.append(self.my_field.e_color)
        if self.my_field.color2 is not None:
            X.append(self.my_field.color2)
            e_X.append(self.my_field.e_color2)
        if self.my_field.plx is not None:
            X.append(self.my_field.plx)
            e_X.append(self.my_field.e_plx)
        if self.my_field.pmra is not None:
            X.append(self.my_field.pmra)
            e_X.append(self.my_field.e_pmra)
        if self.my_field.pmde is not None:
            X.append(self.my_field.pmde)
            e_X.append(self.my_field.e_pmde)
        X = np.array(X)
        e_X = np.array(e_X)

        if X.shape[0] < 3:
            raise ValueError(
                "At least one data dimension beyond (RA, DEC) is required"
                + " for the Bayesian DA method"
            )

        self._vp("{}       : ({:.4f}, {:.4f})".format(cent_str, *xy_center), 1)
        self._vp(f"radius         : {self.my_field.radius:.4f} [deg]", 1)
        self._vp(f"N_cluster      : {self.my_field.N_cluster}", 1)
        self._vp(f"N_runs         : {N_runs}", 1)

        out_mssg, probs = bayesian_mp(
            self.my_field.N_cluster,
            X,
            e_X,
            xy_center,
            self.my_field.radius,
            N_runs,
            self.rng,
        )
        self._vp(out_mssg, 1)
        self.my_field.probs = probs
        self._vp("Probabilities stored in 'probs' attribute", 1)

    def fastmp(
        self,
        fixed_centers: bool = True,
        N_runs: int = 1000,
        eq_to_gal: bool = True,
        seed: int | None = None,
    ) -> None:
        """Assign membership probabilities.

        Estimate the probability of being a true cluster member for all observed
        stars using the fastMP algorithm. The following data dimensions are required:
        ``(pmRA, pmDE, plx)``; photometry is not employed. Center estimates in
        ``(RA, DEC)``, as well as ``(pmRA, pmDE)`` and ``plx`` are required.

        :param fixed_centers: If ``True`` the center values (radec_c, pms_c, plx_c)
            stored in the :py:class:`Cluster <asteca.cluster.Cluster>` object will be
            kept fixed throughout the process
        :type fixed_centers: bool
        :param N_runs: Maximum number of resamples
        :type N_runs: int
        :param eq_to_gal: Convert ``(RA, DEC)`` to ``(lon, lat)``. Useful for clusters
            with large ``DEC`` values to reduce the frame's distortion
        :type eq_to_gal: bool
        :param seed: Random seed. If ``None`` a random integer will be generated and used
        :type seed: int | None

        :raises AttributeError: If the :py:class:`Cluster <asteca.cluster.Cluster>`
            object is missing a required attribute
        :raises ValueError: If the ``N_runs`` parameter is less than 10
        """
        self._vp("\nRunning fastMP")

        # These two are checked separately to avoid type checking issues
        if self.my_field.ra is None or self.my_field.dec is None:
            raise AttributeError("Missing required 'cluster' attributes: (ra, dec)")
        required_attrs = (
            "pmra",
            "pmde",
            "plx",
            "e_pmra",
            "e_pmde",
            "e_plx",
            "radec_c",
            "pms_c",
            "plx_c",
            "N_cluster",
        )
        missing = [
            attr
            for attr in required_attrs
            if not hasattr(self.my_field, attr) or getattr(self.my_field, attr) is None
        ]
        if missing:
            raise AttributeError(
                f"Missing required 'cluster' attributes: {', '.join(missing)}"
            )

        if N_runs < 10:
            raise ValueError("Parameter 'N_runs' should be > 10")

        # Set seed
        self.rng = np.random.default_rng(seed)

        from .modules import cluster_priv as cp
        from .modules.fastmp import fastMP

        xv, yv = self.my_field.ra, self.my_field.dec
        xy_center = self.my_field.radec_c
        cent_str = "radec_c "
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            # Add center as last element
            xv_t = np.array(list(xv) + [xy_center[0]])
            yv_t = np.array(list(yv) + [xy_center[1]])
            # Convert
            xv_t, yv_t = cp.radec2lonlat(xv_t, yv_t)
            # Extract center value
            xv, yv = xv_t[:-1], yv_t[:-1]
            xy_center = (xv_t[-1], yv_t[-1])
            cent_str = "lonlat_c"

        self._vp("{}       : ({:.4f}, {:.4f})".format(cent_str, *xy_center), 1)
        self._vp("pms_c          : ({:.3f}, {:.3f})".format(*self.my_field.pms_c), 1)
        self._vp(f"plx_c          : {self.my_field.plx_c:.3f}", 1)
        self._vp(f"fixed_centers  : {fixed_centers}", 1)
        self._vp(f"N_cluster      : {self.my_field.N_cluster}", 1)
        # centers_ex_flag = False
        # if centers_ex is not None:
        #     centers_ex_flag = True
        # self._vp(f"centers_ex     : {centers_ex_flag}", 1)
        self._vp(f"N_resample     : {N_runs}", 1)

        # Generate input data array for fastMP
        X = np.array(
            [
                xv,
                yv,
                self.my_field.pmra,
                self.my_field.pmde,
                self.my_field.plx,
                self.my_field.e_pmra,
                self.my_field.e_pmde,
                self.my_field.e_plx,
            ]
        )
        out_mssg, probs = fastMP(
            self.rng,
            xy_center,
            self.my_field.pms_c,
            self.my_field.plx_c,
            self.my_field.N_cluster,
            fixed_centers,
            X,
            N_runs,
        )

        self._vp(out_mssg, 1)
        self.my_field.probs = probs
        self._vp("Probabilities stored in 'probs' attribute", 1)

    def select_members(
        self, prob_min: float | None = 0.5, N_prob_max: int | None = None
    ) -> Cluster:
        """Select members based on membership probabilities. This method modifies
        the :py:class:`Cluster <asteca.cluster.Cluster>` object in place, keeping only the
        stars that meet the selection criteria.

        :param prob_min: Minimum probability threshold for membership selection. Stars with
            probabilities above this threshold will be selected as members.
        :type prob_min: float | None
        :param N_prob_max: Maximum number of members to select based on probability. The
            top ``N_prob_max`` stars with the highest probabilities will be selected as
            members.
        :type N_prob_max: int | None

        :raises AttributeError: If the 'probs' attribute is not defined in the
            'cluster' object
        :raises ValueError: Raised if:

            - both ``prob_min`` and ``N_prob_max`` are provided,
            - neither is provided,
            - ``prob_min`` is not between 0 and 1,
            - no members are found with probability above ``prob_min``
            - ``N_prob_max`` is greater than the total number of stars

        :return: A new :py:class:`Cluster <asteca.cluster.Cluster>` object containing
            only the selected members
        :rtype: Cluster
        """
        if self.my_field.probs is None:
            raise AttributeError(
                "The 'probs' attribute is not defined in the 'cluster' object"
            )

        if prob_min is not None and N_prob_max is not None:
            raise ValueError(
                "Only one of 'prob_min' or 'N_prob_max' should be provided"
            )

        N_all = len(self.my_field.probs)

        if N_prob_max is not None:
            if N_prob_max > N_all:
                raise ValueError(
                    f"'N_prob_max' ({N_prob_max}) cannot be greater than the total"
                    + f" number of stars ({N_all})"
                )
            if N_prob_max <= 0:
                raise ValueError("'N_prob_max' must be a positive integer")

            # Sort probabilities and select top N_prob_max
            sorted_indices = np.argsort(self.my_field.probs)[::-1]
            msk = np.zeros_like(self.my_field.probs, dtype=bool)
            msk[sorted_indices[:N_prob_max]] = True
            self._vp(f"\nTop {N_prob_max} members selected based on probability")
        elif prob_min is not None:
            if not (0.0 <= prob_min <= 1.0):
                raise ValueError("'prob_min' should be between 0 and 1")
            msk = self.my_field.probs >= prob_min
            if not msk.any():
                raise ValueError(f"No members found with probability >= {prob_min}")
            self._vp(f"\nMembers selected with probability >= {prob_min}")
        else:
            raise ValueError("Either 'prob_min' or 'N_prob_max' must be provided")

        from .cluster import Cluster

        new_cluster = Cluster(verbose=0)

        # Apply msk to self object
        for attr, arr in vars(self.my_field).items():
            # Only apply mask to array attributes with the same length as N_stars
            if isinstance(arr, np.ndarray) and len(arr) == self.my_field.N_stars:
                setattr(new_cluster, attr, arr[msk])
            else:
                setattr(new_cluster, attr, arr)

        new_cluster.N_stars = msk.sum()
        new_cluster.N_cluster = new_cluster.N_stars
        self._vp("New Cluster object created", 1)
        self._vp(f"N_cluster      : {new_cluster.N_cluster} (kept)", 1)
        self._vp(f"Field stars    : {N_all - new_cluster.N_cluster} (removed)", 1)

        return new_cluster
