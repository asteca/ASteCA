import os
import numpy as np
from dataclasses import dataclass
from .modules import isochrones_priv


@dataclass
class Isochrones:
    """Define an :py:class:`Isochrones` object.

    This object contains the loaded theoretical isochrones used by the
    :py:class:`Synthetic <asteca.synthetic.Synthetic>` class to generate synthetic
    clusters. See :ref:`isochronesload` for more details.

    :param model: The model must be one of the supported isochrone services:
        `PARSEC <http://stev.oapd.inaf.it/cgi-bin/cmd_3.7>`__,
        `MIST <https://waps.cfa.harvard.edu/MIST/>`__, or
        `BASTI <http://basti-iac.oa-abruzzo.inaf.it/isocs.html>`__.
    :type model: str: ``PARSEC, MIST, BASTI``
    :param isochs_path: Path to the folder that contains the files for the theoretical
        isochrones
    :type isochs_path: str
    :param magnitude: Magnitude's filter name as defined in the theoretical isochrones.
        Example for Gaia's ``G`` magnitude: ``"Gmag"``
    :type magnitude: str
    :param color: Tuple containing the filter names that generate the first color
        defined. The correct format is: ``("filter1", "filter2")``, where ``filter1``
        and ``filter2`` are the names of the filters that are combined to generate
        the color. The order is important because the color will be generated as:
        ``filter1-filter2``. Example for Gaia's 'BP-RP' color:
        ``("G_BPmag", "G_RPmag")``
    :type color: tuple
    :param color2: Second color to use in the analysis. Same format as that
        used by the ``color`` parameter, defaults to ``None``
    :type color2: tuple | None
    :param magnitude_effl: Effective lambda (in Angstrom) for the magnitude filter,
        defaults to ``None``
    :type magnitude_effl: float | None
    :param color_effl: Effective lambdas for the filters that make up the ``color``
        defined in the :py:class:`Isochrones` object. E.g.:
        ``(1111.11, 2222.22)`` where ``1111.11`` and ``2222.22`` are the effective
        lambdas (in Angstrom) for each filter, in the same order as ``color``, defaults
        to ``None``
    :type color_effl: tuple | None
    :param color2_effl: Same as ``color_effl`` but for a second (optional) color
        defined, defaults to ``None``
    :type color2_effl: tuple | None
    :param z_to_FeH: If ``None``, the default ``z`` values in the isochrones will be
        used to generate the synthetic clusters. If ``float``, it must represent the
        solar metallicity for these isochrones. The metallicity values will then be
        converted to ``[FeH]`` values, to be used by the
        :py:meth:`Synthetic.generate() <asteca.synthetic.Synthetic.generate>` method,
        defaults to ``None``
    :type z_to_FeH: float | None
    :param column_names: Column names for the initial mass, metallicity, and age for
        the photometric system's isochrones files. Example:
        ``{"mass_col": "Mini", "met_col": "Zini", "age_col": "logAge"}``.
        This dictionary is defined internally in **ASteCA** and should only be given
        by the user if the isochrone service changes its format and the `isochrones`
        class fails to load the files, defaults to ``None``
    :type column_names: dict, optional
    :param N_interp: Number of interpolation points used to ensure that all isochrones
        are the same shape, defaults to ``2500``
    :type N_interp: int
    :param parsec_rm_stage_9: If the isochrones are PARSEC, this argument set to
        ``True`` will remove the *post_AGB* stage (label=9) which are still
        "`in preparation <http://stev.oapd.inaf.it/cmd_3.7/faq.html>`__", defaults
        to ``True``
    :type parsec_rm_stage_9: bool
    """

    model: str
    isochs_path: str
    magnitude: str
    color: tuple
    color2: tuple | None = None
    magnitude_effl: float | None = None
    color_effl: tuple | None = None
    color2_effl: tuple | None = None
    z_to_FeH: float | None = None
    column_names: dict | None = None
    N_interp: int = 2500
    parsec_rm_stage_9: bool = True

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

        print("\nInstantiating isochrones...")
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

        N_met, N_age, _, N_isoch = self.theor_tracks.shape
        print(f"Model          : {self.model}")
        print(f"N_files        : {N_isoch_files}")
        print(f"N_met          : {N_met}")
        print(f"N_age          : {N_age}")
        print(f"N_isoch        : {N_isoch}")
        print(f"{met_n} range      : [{zmin}, {zmax}]")
        print(f"loga range     : [{amin}, {amax}]")
        print(f"Magnitude      : {self.magnitude}")
        print(f"Color          : {self.color[0]}-{self.color[1]}")
        if self.color2 is not None:
            print(f"Color2         : {self.color2[0]}-{self.color2[1]}")
        print("Isochrone object generated")

    def _func_z_to_FeH(self, z_to_FeH):
        """Convert z to FeH"""
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
        """Return the minimum and maximum values for the metallicity and age defined
        in the theoretical isochrones.

        :return: Tuple of (minimum_metallicity, maximum_metallicity, minimum_age,
            maximum_age)
        :rtype: tuple[float]
        """
        zmin = self.met_age_dict["met"].min()
        zmax = self.met_age_dict["met"].max()
        amin = self.met_age_dict["loga"].min()
        amax = self.met_age_dict["loga"].max()
        return zmin, zmax, amin, amax
