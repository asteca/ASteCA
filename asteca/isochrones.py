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
    See :ref:`isochronesload` for more details.

    Parameters
    ----------
    model : str, {"PARSEC", "MIST", "BASTI"}
        The model must be one of the supported isochrone services:
        `PARSEC <http://stev.oapd.inaf.it/cgi-bin/cmd_3.7>`_,
        `MIST <https://waps.cfa.harvard.edu/MIST/>`_, or
        `BASTI <http://basti-iac.oa-abruzzo.inaf.it/isocs.html>`_.
    isochs_path : str
        Path to the folder that contains the files for the theoretical isochrones.
    magnitude : str
        Magnitude's filter name as defined in the theoretical isochrones.
        Example for Gaia's ``G`` magnitude: ``"Gmag"``.
    color : tuple
        Tuple containing the filter names that generate the first color defined.
        The correct format is: ``("filter1", "filter2")``, where ``filter1``
        and ``filter2`` are the names of the filters that are combined to generate
        the color. The order is important because the color will be generated as:
        ``filter1-filter2``. Example for Gaia's 'BP-RP' color:
        ``("G_BPmag", "G_RPmag")``.
    color2 : tuple, optional, default=None
        Optional second color to use in the analysis. Same format as that used by the
        ``color`` parameter.
    magnitude_effl : float, optional, default=None
        Effective lambda (in Angstrom) for the magnitude filter.
    color_effl : tuple, optional, default=None
        Effective lambdas for the filters that make up the ``color`` defined in the
        :py:mod:`asteca.isochrones` object. E.g.: ``(1111.11, 2222.22)`` where
        ``1111.11`` and ``2222.22`` are the effective lambdas (in Angstrom) for each
        filter, in the same order as ``color``.
    color2_effl : tuple, optional, default=None
        Same as ``color_effl`` but for a second (optional) color defined.
    z_to_FeH : float, optional, default=None
        If ``None``, the default ``z`` values in the isochrones will be used to
        generate the synthetic clusters. If ``float``, it must represent the solar
        metallicity for these isochrones. The metallicity values will then be converted
        to ``[FeH]`` values, to be used by the :meth:`synthetic.generate()` method.
    N_interp : int, default=2500
        Number of interpolation points used to ensure that all isochrones are the
        same shape.
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

    model: str
    isochs_path: str
    magnitude: str
    color: tuple
    color2: Optional[tuple] = None
    magnitude_effl: Optional[float] = None
    color_effl: Optional[tuple] = None
    color2_effl: Optional[tuple] = None
    N_interp: int = 2500
    z_to_FeH: Optional[float] = None
    column_names: Optional[dict] = None
    parsec_rm_stage_9: Optional[bool] = True

    def __post_init__(self):
        # Check that the number of colors match
        if self.color2 is not None and self.color2_effl is None:
            raise ValueError(
                "Second color is defined but its effective lambdas are missing."
            )
        if self.color2 is None and self.color2_effl is not None:
            raise ValueError(
                "Lambdas for the second color are defined but second color is missing."
            )

        # Check model input
        self.model = self.model.upper()
        models = ("PARSEC", "MIST", "BASTI")
        if self.model not in models:
            raise ValueError(
                f"Model '{self.model}' not recognized. Should be one of {models}"
            )

        # Check path to isochrones
        if os.path.isdir(self.isochs_path) is False:
            raise ValueError(f"Path '{self.isochs_path}' not found")

        print("Instantiating isochrones...")
        # Load isochrone files
        self.theor_tracks, self.color_filters, self.met_age_dict, N_isoch_files = (
            isochrones_priv.load(self)
        )

        # Convert z to FeH if requested
        met_n = "z  "
        if self.z_to_FeH is not None:
            self._func_z_to_FeH(self.z_to_FeH)
            met_n = "FeH"
        zmin, zmax, amin, amax = self._min_max()

        N_met, N_age, N_cols, N_interp = self.theor_tracks.shape
        print(f"Model          : {self.model}")
        print(f"N_files        : {N_isoch_files}")
        print(f"N_met          : {N_met}")
        print(f"N_age          : {N_age}")
        print(f"N_interp       : {N_interp}")
        print(f"{met_n} range      : [{zmin}, {zmax}]")
        print(f"loga range     : [{amin}, {amax}]")
        print(f"Magnitude      : {self.magnitude}")
        print(f"Color          : {self.color[0]}-{self.color[1]}")
        if self.color2 is not None:
            print(f"Color2         : {self.color2[0]}-{self.color2[1]}")
        print("Isochrone object generated\n")

    def _func_z_to_FeH(self, z_to_FeH):
        r"""Convert z to FeH"""
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
