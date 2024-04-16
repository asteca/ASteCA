import os
from dataclasses import dataclass
from typing import Optional
from .modules import isochrones_priv


@dataclass
class isochrones:
    r"""Define an ``isochrones`` object.

    This object contains the loaded theoretical isochrones used by the
    :class:`synthetic` class to generate synthetic clusters.

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
        print("Isochrone object generated\n")
