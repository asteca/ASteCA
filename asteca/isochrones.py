import numpy as np


class Isochrones:
    """Define an :py:class:`Isochrones` object.

    This object contains the loaded theoretical isochrones used by the
    :py:class:`Synthetic <asteca.synthetic.Synthetic>` class to generate synthetic
    clusters. See :ref:`isochrones_module` for more details.

    :param model: The model must be one of the supported isochrone services:
        `PARSEC <http://stev.oapd.inaf.it/cgi-bin/cmd_3.7>`__,
        `MIST <https://waps.cfa.harvard.edu/MIST/>`__, or
        `BASTI <http://basti-iac.oa-abruzzo.inaf.it/isocs.html>`__.
    :type model: str
    :param isochs_path: Path to the file or folder that contains the files for the
        theoretical isochrones
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
        :py:meth:`Synthetic.generate` method, defaults to ``None``
    :type z_to_FeH: float | None
    :param N_interp: Number of interpolation points used to ensure that all isochrones
        are the same shape, defaults to ``2000``
    :type N_interp: int
    :param parsec_rm_stage_9: If the isochrones are PARSEC, this argument set to
        ``True`` will remove the *post_AGB* stage (label=9) which are still
        "`in preparation <http://stev.oapd.inaf.it/cmd/faq.html>`__", defaults
        to ``True``
    :type parsec_rm_stage_9: bool
    :param column_names: Column names for the initial mass, metallicity, and age for
        the photometric system's isochrones files. Example:
        ``{"mass_col": "Mini", "met_col": "Zini", "age_col": "logAge"}``.
        This dictionary is defined internally in **ASteCA** and should only be given
        by the user if the isochrone service changes its format and the `isochrones`
        class fails to load the files, defaults to ``None``
    :type column_names: dict | None
    :param verbose: Verbose level. A value of ``0`` hides all output, defaults to ``1``
    :type verbose: int

    :raises ValueError: If any of the attributes is not recognized as a valid option,
        or there are missing required attributes
    """

    def __init__(
        self,
        model: str,
        isochs_path: str,
        magnitude: str,
        color: tuple,
        color2: tuple | None = None,
        magnitude_effl: float | None = None,
        color_effl: tuple | None = None,
        color2_effl: tuple | None = None,
        z_to_FeH: float | None = None,
        N_interp: int = 2000,
        parsec_rm_stage_9: bool = True,
        column_names: dict | None = None,
        verbose: int = 1,
    ) -> None:
        self.model = model
        self.isochs_path = isochs_path
        self.magnitude = magnitude
        self.color = color
        self.color2 = color2
        self.magnitude_effl = magnitude_effl
        self.color_effl = color_effl
        self.color2_effl = color2_effl
        self.z_to_FeH = z_to_FeH
        self.column_names = column_names
        self.N_interp = N_interp
        self.parsec_rm_stage_9 = parsec_rm_stage_9
        self.verbose = verbose

        from .modules import isochrones_priv

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

        self._vp("\nInstantiating isochrones")
        # Load isochrone files
        self.theor_tracks, self.color_filters, self.met_age_dict, N_isoch_files = (
            isochrones_priv.load(
                self.model,
                self.isochs_path,
                self.magnitude,
                self.color,
                self.color2,
                self.column_names,
                self.N_interp,
                self.parsec_rm_stage_9,
            )
        )

        # Convert z to FeH if requested
        met_n = "z  "
        if self.z_to_FeH is not None:
            self._func_z_to_FeH(self.z_to_FeH)
            met_n = "FeH"

        # Extract metallicity and age ranges
        self.zmin = self.met_age_dict["met"].min()
        self.zmax = self.met_age_dict["met"].max()
        self.amin = self.met_age_dict["loga"].min()
        self.amax = self.met_age_dict["loga"].max()

        N_met, N_age, _, N_isoch = self.theor_tracks.shape
        self._vp(f"Model          : {self.model}", 1)
        self._vp(f"N_files        : {N_isoch_files}", 1)
        self._vp(f"N_mets         : {N_met}", 1)
        self._vp(f"N_ages         : {N_age}", 1)
        self._vp(f"N_isochs       : {N_isoch}", 1)
        self._vp(f"{met_n}  range     : [{self.zmin}, {self.zmax}]", 1)
        self._vp(f"loga range     : [{self.amin}, {self.amax}]", 1)
        self._vp(f"Magnitude      : {self.magnitude}", 1)
        self._vp(f"Color          : {self.color[0]}-{self.color[1]}", 1)
        if self.color2 is not None:
            self._vp(f"Color2         : {self.color2[0]}-{self.color2[1]}", 1)
        self._vp("Isochrone object generated")

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

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

    def _interpolate(self, met: float, loga: float) -> np.ndarray:
        """ """
        met_msk = self.met_age_dict["met"] == met
        age_msk = self.met_age_dict["loga"] == loga

        # Index of the first element in mask that is True
        try:
            met_i = np.where(met_msk)[0][0]
            age_i = np.where(age_msk)[0][0]
            isoch = self.theor_tracks[met_i][age_i]
        except IndexError:
            from .modules import synth_cluster_priv as scp

            params = {
                "met": met,
                "loga": loga,
                "alpha": 0,
                "beta": 0,
                "Av": 0,
                "DR": 0,
                "Rv": 0,
                "dm": 0,
            }
            ml, mh, al, ah = scp.properModel(self.met_age_dict, {}, params)[-4:]
            # Add dimension of zeroes to array
            padded_arr = np.pad(self.theor_tracks, ((0, 0), (0, 0), (0, 2), (0, 0)))

            m_ini_idx = self.theor_tracks.shape[2]

            isoch = scp.zaWAverage(
                padded_arr,
                self.met_age_dict,
                m_ini_idx,
                met,
                loga,
                ml,
                mh,
                al,
                ah,
            )[:m_ini_idx, :]

        return isoch
