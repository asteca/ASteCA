import os
import numpy as np
from dataclasses import dataclass
from typing import Optional
from .modules import isochrones_priv


@dataclass
class isochrones:
    r"""Define an ``isochrones`` object.

    This object contains the loaded theoretical isochrones used by the
    :py:mod:`asteca.synthetic` class to generate synthetic clusters.

    Parameters
    ----------
    isochs_path : str
        Path to the folder that contains the files for the theoretical isochrones.
        The name of the folder must be one of the supported isochrone services:
        `PARSEC <http://stev.oapd.inaf.it/cgi-bin/cmd_3.7>`_,
        `MIST <https://waps.cfa.harvard.edu/MIST/>`_, or
        `BASTI <http://basti-iac.oa-abruzzo.inaf.it/isocs.html>`_.
        Examples of valid paths: ``isochrones/PARSEC/``, ``Mist/``, ``basti``.
        See :ref:`isochronesload` for more detailed information on how to properly
        store the isochrone files.
    magnitude : dict
        Dictionary containing the magnitude's filter name (as defined in the files
        of the theoretical isochrones) as the key, and its effective lambda
        (in Angstrom) as the value. Example for Gaia's ``G`` magnitude:
        ``{"Gmag": 6390.21}``.
    color : dict
        Dictionary containing the color used in the cluster's analysis. The correct
        format is: ``{"filter1": 1111.11, "filter2": 2222.22}``, where ``filter1``
        and `filter2` are the names of the filters that are combined to generate
        the color. The order is important because the color will be generated as:
        ``filter1-filter2``. The values ``1111.11`` and ``2222.22`` are the effective
        lambdas (in Angstrom) for each filter. The color does not need to be defined
        in the same photometric system as the magnitude.
        Example for Gaia's 'BP-RP' color:
        ``{"G_BPmag": 5182.58, "G_RPmag": 7825.08}``
    z_to_FeH : float, optional, default=None
        If ``None``, the default ``z`` values (defined when loading the isochrones
        via the :py:mod:`asteca.isochrones` object) will be used to generate the
        synthetic clusters. If ``float``, it must represent the solar metallicity
        for these isochrones. The metallicity values will then be converted to
        ``[FeH]`` values, to be used by the :meth:`synthetic.generate()` method.
    N_interp : int, default=2500
        Number of interpolation points used to ensure that all isochrones are the
        same shape.
    color2 : dict, optional, default=None
        Optional second color to use in the analysis. Same format as that used by the
        ``color`` parameter.
    column_names : dict, optional, default=None
        Column names for the initial mass, metallicity, and age for the photometric
        system's isochrones files. Example:
        ``{"mass_col": "Mini", "met_col": "Zini", "age_col": "logAge"}``.
        This dictionary is defined internally in `ASteCA` and should only be given
        by the user if the isochrone service changes its format and the `isochrones`
        class fails to load the files.
    parsec_rm_stage_9 : boll, optional, default=True
        If the isochrones are PARSEC, this argument set to ``True`` will remove the
        *post_AGB* stage (label=9) which are still
        "`in preparation <http://stev.oapd.inaf.it/cmd_3.7/faq.html>`_".

    """
    isochs_path: str
    magnitude: dict
    color: dict
    z_to_FeH: float | None = None
    N_interp: int = 2500
    color2: Optional[dict] = None
    column_names: Optional[dict] = None
    parsec_rm_stage_9: Optional[bool] = True

    def __post_init__(self):
        # Extract model from isochrones' path
        if self.isochs_path.endswith("/"):
            self.model = self.isochs_path[:-1].split("/")[-1].upper()
        else:
            self.model = self.isochs_path.split("/")[-1].upper()

        # Check model input
        models = ("PARSEC", "MIST", "BASTI")
        if self.model not in models:
            raise ValueError(
                f"Model '{self.model}' not recognized. Should be one of {models}"
            )

        # Check path to isochrones
        if os.path.isdir(self.isochs_path) is False:
            raise ValueError(
                f"Path '{self.isochs_path}' must point to the folder that "
                + f"contains the '{self.model}' isochrones"
            )

        # Extract magnitude, color(s), and lambdas
        self.mag_color_lambdas = self.magnitude | self.color
        self.mag_filter_name = list(self.magnitude.keys())[0]
        self.color_filter_name = [list(self.color.keys())]
        # Add second color if any
        if self.color2 is not None:
            self.color_filter_name.append(list(self.color2.keys()))
            self.mag_color_lambdas = self.mag_color_lambdas | self.color2

        # Load isochrone files
        self.theor_tracks, self.color_filters, self.met_age_dict = isochrones_priv.load(
            self
        )

        # Convert z to FeH if requested
        met_n = 'z  '
        if self.z_to_FeH is not None:
            self._func_z_to_FeH(self.z_to_FeH)
            met_n = 'FeH'

        zmin, zmax, amin, amax = self._min_max()
        print(f"{met_n} range  : [{zmin}, {zmax}]")
        print(f"loga range : [{amin}, {amax}]")
        print("Isochrone object generated\n")

    def _func_z_to_FeH(self, z_to_FeH):
        r"""Convert z to FeH
        """
        feh = np.log10(self.met_age_dict["met"] / z_to_FeH)
        N_old = len(feh)
        round_n = 4
        while True:
            feh_r = np.round(feh, round_n)
            N_new = len(set(feh_r))
            # If no duplicated values exist after rounding
            if N_old == N_new:
                break
            round_n += 1
        # Replace old values
        self.met_age_dict["met"] = feh_r

    def _min_max(self) -> tuple[float]:
        r"""Return the minimum and maximum values for the metallicity and age defined
        in the theoretical isochrones.

        Returns
        -------
        tuple[float]
            Tuple of (minimum_metallicity, maximum_metallicity, minimum_age, maximum_age)

        """
        zmin = self.met_age_dict["met"].min()
        zmax = self.met_age_dict["met"].max()
        amin = self.met_age_dict["loga"].min()
        amax = self.met_age_dict["loga"].max()
        return zmin, zmax, amin, amax
