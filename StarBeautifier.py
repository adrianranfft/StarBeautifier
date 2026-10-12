"""
StarBeautifier  -  Siril-Script (sirilpy, Siril >= 1.4)
(vormals "StarBeautyfier" / "Star Protect Reduce")
=============================================================

Sternreduktion mit Katalog-Schutzmaske: Sterne heller als eine
Magnitudengrenze (Hipparcos-Katalog, online über VizieR) bleiben
unverändert, alle anderen Sterne werden gedämpft.

Voraussetzung:
  - Geladenes Bild = Sternbild aus StarNet (oder jedes andere Sternbild)
  - Das Bild ist plate-gesolvt (Werkzeuge -> Astrometrie -> Plate Solving,
    möglichst mit Verzeichnungskorrektur / SIP). Ein reines Sternbild lässt
    sich problemlos solven. Alternativ im Dialog "WCS aus Datei" wählen.

Ergebnis:
  - Das geladene Bild wird durch das reduzierte Sternbild ersetzt (Undo möglich)
  - Neben dem Bild werden gespeichert:
      <name>_protect_mask.fit   Schutzmaske 0..1
      <name>_protect_check.jpg  Kontrollbild (grün = geschützt, rot = nicht gefunden)
      <name>_protect_stars.csv  Katalogsterne mit Position und Radius

Rechnung:
  L_neu    = L ^ Exponent                  (farberhaltend über Luminanz)
  ergebnis = maske * stars + (1 - maske) * faktor * stars_reduziert

Danach z. B. in PixelMath:  starless + stars
"""

import os
import sys
import csv
import json

import sirilpy as s
from sirilpy import tksiril, SirilError

s.ensure_installed("ttkthemes")
s.ensure_installed("scipy")
s.ensure_installed("astropy")
s.ensure_installed("astroquery")
s.ensure_installed("pillow")

import numpy as np  # noqa: E402
import tkinter as tk  # noqa: E402
from tkinter import ttk, filedialog, messagebox  # noqa: E402
from ttkthemes import ThemedTk  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
import warnings  # noqa: E402
from astropy.wcs import FITSFixedWarning  # noqa: E402
from astropy.utils.exceptions import AstropyWarning  # noqa: E402

# astropy-Hinweise (MJD-OBS ergänzt, Header ohne NAXIS) sind harmlos -> nicht ins Siril-Log
warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", AstropyWarning)

# =============================================================================
# Sprache / Language
# =============================================================================

LANG = "de"

EN = {
    # --- Kern / Log
    "Feldmitte RA {ra:.2f}°, Dec {dec:.2f}°, Radius {rad:.1f}°":
        "Field centre RA {ra:.2f}°, Dec {dec:.2f}°, radius {rad:.1f}°",
    "Frage Hipparcos ab (V < {m}) ...": "Querying Hipparcos (V < {m}) ...",
    "{n} Katalogsterne heller als {m} mag im Bild": "{n} catalogue stars brighter than {m} mag in the image",
    "Hinweis: Hipparcos ist nur bis ca. 7.5–8 mag vollständig.":
        "Note: Hipparcos is only complete to about 7.5–8 mag.",
    "Keine Katalogsterne im Bildfeld gefunden – Plate-Solve prüfen.":
        "No catalogue stars found in the field – check the plate solve.",
    "Orientierung: {o} (enge Treffer {a} / {b} von {n}, Median {ma:.1f} / {mb:.1f} px)":
        "Orientation: {o} (close matches {a} / {b} of {n}, median {ma:.1f} / {mb:.1f} px)",
    "Gefunden: {f}, nicht gefunden: {m}, Positionsabweichung Median {d:.1f} px":
        "Found: {f}, not found: {m}, median position offset {d:.1f} px",
    "WARNUNG: große Positionsabweichung – Plate-Solve (Verzeichnung/SIP) prüfen oder Suchradius erhöhen.":
        "WARNING: large position offset – check the plate solve (distortion/SIP) or increase the search radius.",
    "Für 'Glow nur für Katalogsterne' wird die Sternliste eines vorherigen Laufs benötigt – zuerst 'Anwenden' oder 'Vorschau'.":
        "'Glow catalogue stars only' needs the star list of a previous run – use 'Apply' or 'Preview' first.",
    "Gespeichert: {b}_protect_mask.fit / _check.jpg / _stars.csv":
        "Saved: {b}_protect_mask.fit / _check.jpg / _stars.csv",
    # --- GUI
    "Fehler": "Error",
    "Keine Verbindung zu Siril.": "No connection to Siril.",
    "Kein Bild geladen. Bitte zuerst das Sternbild öffnen.": "No image loaded. Please open the star image first.",
    "Helle Katalogsterne schützen, Rest reduzieren": "Protect bright catalogue stars, reduce the rest",
    "Glow nur für Katalogsterne (Mag-Bereich)": "Glow catalogue stars only (mag range)",
    "Nur Glow anwenden": "Apply glow only",
    "Maske + Kontrollbild speichern": "Save mask + check image",
    "WCS: aus geladenem Bild": "WCS: from loaded image",
    "WCS aus Datei…": "WCS from file…",
    "Nur nötig, wenn das geladene Bild nicht plate-gesolvt ist: geplatesolvtes FITS gleicher Größe wählen.":
        "Only needed if the loaded image is not plate-solved: choose a plate-solved FITS of the same size.",
    "Voreinstellung": "Preset",
    "Preset wählen oder neuen Namen eintippen.": "Choose a preset or type a new name.",
    "Speichern": "Save",
    "Löschen": "Delete",
    "Werkseinstellungen": "Factory defaults",
    "Brennweite laut Bild: unbekannt": "Focal length from image: unknown",
    "Brennweite laut Bild: {f:.0f} mm": "Focal length from image: {f:.0f} mm",
    "Schließen": "Close",
    "Vorschau…": "Preview…",
    "Anwenden": "Apply",
    "Sprache": "Language",
    "Die gewählte Datei enthält kein WCS.": "The selected file contains no WCS.",
    "Größe passt nicht: WCS-Datei {a}x{b}, Bild {c}x{d}": "Size mismatch: WCS file {a}x{b}, image {c}x{d}",
    "Das Bild ist nicht plate-gesolvt. Bitte zuerst Plate Solving ausführen oder 'WCS aus Datei…' wählen.":
        "The image is not plate-solved. Please run plate solving first or choose 'WCS from file…'.",
    "StarBeautifier: Brennweite aus Plate-Solve geschätzt: {f:.1f} mm":
        "StarBeautifier: focal length estimated from plate solve: {f:.1f} mm",
    "StarBeautifier: WCS gelesen, SIP-Verzeichnung: {s}": "StarBeautifier: WCS read, SIP distortion: {s}",
    "ja": "yes",
    "nein": "no",
    "Glow-Stärke ist 0 – bitte einen Wert > 0 einstellen.": "Glow strength is 0 – please set a value > 0.",
    "StarBeautifier: Glow {g} (Größe {s}, Weichheit {w}) nachträglich angewendet":
        "StarBeautifier: glow {g} (size {s}, softness {w}) applied afterwards",
    "StarBeautifier: Glow {g} (Größe {s}, Weichheit {w}) angewendet":
        "StarBeautifier: glow {g} (size {s}, softness {w}) applied",
    "Erklärung": "Explanation",
    "Übersicht": "Overview",
    "StarBeautifier: Brennweite aus Bild: {f:.1f} mm": "StarBeautifier: focal length from image: {f:.1f} mm",
    "Kein Preset namens „{n}“.": "No preset named “{n}”.",
    "StarBeautifier: Preset „{n}“ geladen": "StarBeautifier: preset “{n}” loaded",
    "Bitte einen Namen für das Preset eintippen.": "Please type a name for the preset.",
    "Preset „{n}“ überschreiben?": "Overwrite preset “{n}”?",
    "StarBeautifier: Preset „{n}“ gespeichert ({p})": "StarBeautifier: preset “{n}” saved ({p})",
    "Preset „{n}“ löschen?": "Delete preset “{n}”?",
    "StarBeautifier: Preset „{n}“ gelöscht": "StarBeautifier: preset “{n}” deleted",
    "StarBeautifier: Werkseinstellungen geladen": "StarBeautifier: factory defaults loaded",
    "StarBeautifier: fertig. Erneutes 'Anwenden' rechnet wieder vom Original.":
        "StarBeautifier: done. Pressing 'Apply' again starts from the original.",
    "StarBeautifier Fehler: {e}": "StarBeautifier error: {e}",
    "StarBeautifier: lese Bild …": "StarBeautifier: reading image …",
    "StarBeautifier: Katalog & Sternsuche …": "StarBeautifier: catalogue & star search …",
    "StarBeautifier: Reduktion …": "StarBeautifier: reduction …",
    "StarBeautifier: Glow …": "StarBeautifier: glow …",
    "StarBeautifier: schreibe Ergebnis …": "StarBeautifier: writing result …",
    # --- Vorschau
    "StarBeautifier – Vorschau": "StarBeautifier – Preview",
    "Klick ins Bild wählt\nden Ausschnitt.": "Click into the image\nto choose the crop.",
    "Aktualisieren": "Refresh",
    "Vorher": "Before",
    "Nachher": "After",
    "rechne …": "computing …",
    "StarBeautifier Vorschau-Fehler: {e}": "StarBeautifier preview error: {e}",
    "1:1": "1:1",
    "2:1 (Zoom)": "2:1 (zoom)",
    "1:2 (Übersicht)": "1:2 (overview)",
    # --- Parameter (Label + Tooltip)
    "Schutz bis Magnitude": "Protect up to magnitude",
    "Sterne heller als dieser V-Wert bleiben unverändert (Hipparcos, vollständig bis ca. 7.5).":
        "Stars brighter than this V magnitude stay unchanged (Hipparcos, complete to about 7.5).",
    "Exponent (Stärke)": "Exponent (strength)",
    "Dämpfung der übrigen Sterne: L^Exponent. 2 = moderat, 3 = stark.":
        "Damping of the remaining stars: L^exponent. 2 = moderate, 3 = strong.",
    "Helligkeitsfaktor": "Brightness factor",
    "Zusätzlicher Faktor für die übrigen Sterne.": "Additional factor for the remaining stars.",
    "Maskenrand (px)": "Mask margin (px)",
    "Vergrößert die Schutzkreise, falls Halos angeschnitten werden.":
        "Enlarges the protection circles if halos get clipped.",
    "Weiche Kante (px)": "Feather (px)",
    "Breite des weichen Übergangs am Maskenrand.": "Width of the soft transition at the mask edge.",
    "Suchradius (px)": "Search radius (px)",
    "Wie weit vom Katalogort nach dem Stern gesucht wird (Restverzeichnung).":
        "How far from the catalogue position to search for the star (residual distortion).",
    "Randschwelle": "Edge threshold",
    "Helligkeit, ab der der Sternrand endet (für die Radiusmessung).":
        "Brightness at which the star edge ends (for radius measurement).",
    "Halo-Übergang": "Halo transition",
    "Schutz folgt dem Sternprofil: unterhalb Randschwelle x Wert wird auch im Schutzkreis weich reduziert. Höher = weicherer Übergang, 1 = harter Kreis (alte Version).":
        "Protection follows the star profile: below edge threshold × value it is softly reduced. Higher = softer transition, 1 = hard circle (old version).",
    "Glow-Stärke": "Glow strength",
    "Zusätzlicher weicher Schein um die Sterne. 0 = aus.": "Additional soft glow around the stars. 0 = off.",
    "Glow-Größe": "Glow size",
    "Größe des Scheins relativ zum Sternradius.": "Size of the glow relative to the star radius.",
    "Glow-Weichheit": "Glow softness",
    "0 = kompakter Schein, 1 = langer, weicher Auslauf.": "0 = compact glow, 1 = long, soft falloff.",
    "Glow von Mag": "Glow from mag",
    "Hellste Sterne, die einen Glow bekommen (kleinere Zahl = heller).":
        "Brightest stars that get a glow (smaller number = brighter).",
    "Glow bis Mag": "Glow to mag",
    "Schwächste Sterne, die einen Glow bekommen.": "Faintest stars that get a glow.",
    "Glow-Abnahme": "Glow falloff",
    "0 = alle gleich, 1 = linear von 'Glow von' (voll) bis 'Glow bis' (null), >1 = schneller.":
        "0 = all equal, 1 = linear from 'Glow from' (full) to 'Glow to' (zero), >1 = faster.",
}


def _(text, **kw):
    """Übersetzt einen (deutschen) Text in die aktive Sprache und formatiert ihn."""
    if LANG == "en":
        text = EN.get(text, text)
    return text.format(**kw) if kw else text


def detect_language():
    try:
        import locale
        loc = (locale.getlocale()[0] or "") + (os.environ.get("LANG") or "")
        if "de" in loc.lower() or "german" in loc.lower():
            return "de"
    except Exception:
        pass
    return "en"


HELP_OVERVIEW_EN = (
    "StarBeautifier\n\n"
    "Bright stars from the Hipparcos catalogue are protected, all other stars in the "
    "star image are dimmed. Constellations stay recognisable while the dense star field "
    "recedes and nebulae and dark clouds get more room.\n\n"
    "Workflow:\n"
    "1. Load the star image from StarNet (not the combined image).\n"
    "2. Plate solve with distortion correction (SIP).\n"
    "3. Adjust the values, Apply, check the check image.\n"
    "4. Change values and apply again – it always starts from the original.\n"
    "5. Then in PixelMath: starless + stars\n\n"
    "Hover over a parameter to see its explanation here."
)

HELP_EN = {
    "maglim": (
        "Protect up to magnitude",
        "Stars brighter than this value (V magnitude from Hipparcos) stay unchanged. "
        "Smaller number = brighter star.\n\n"
        "Guide values:\n"
        "• 4–5: only the prominent constellation stars\n"
        "• 5.5–6: everything visible to the naked eye under a dark sky\n"
        "• 7–8: very many protected stars – the reduction barely works and the field "
        "can look uneven\n\n"
        "Focal length: for wide angle (14–24 mm) rather 5–5.5, otherwise very many stars "
        "are protected. For normal and tele lenses (50 mm+) the field is smaller, 6–7 fits.\n\n"
        "Hipparcos is complete to about 7.5 mag."
    ),
    "power": (
        "Exponent (strength)",
        "How strongly the unprotected stars are dimmed. The brightness is raised to the "
        "power of the exponent – faint stars lose disproportionately, bright ones only "
        "a little. The star colour is preserved.\n\n"
        "Example: pixel brightness 0.3\n"
        "• exponent 1.5 → 0.16\n"
        "• exponent 2 → 0.09\n"
        "• exponent 3 → 0.03\n\n"
        "Guide values: 1.5 subtle, 2 moderate, 2.5–3 strong. Above 3 the field quickly "
        "looks empty and unnatural."
    ),
    "dim": (
        "Brightness factor",
        "Additional linear factor for all unprotected stars, applied after the exponent. "
        "1.0 = no additional dimming.\n\n"
        "Difference to the exponent: the factor dims all stars evenly, the exponent mainly "
        "the faint ones. With 0.7–0.8 the medium-bright stars can be toned down too without "
        "raising the exponent further."
    ),
    "margin": (
        "Mask margin (px)",
        "Enlarges the measured protection circle of each bright star by this many pixels.\n\n"
        "Increase if the edges of bright stars look clipped or frayed. Since protection "
        "follows the star profile, a large margin is uncritical.\n\n"
        "Focal length: with longer focal lengths and higher resolution stars get larger in "
        "pixels – then rather 5–10 instead of 3."
    ),
    "feather": (
        "Feather (px)",
        "Width of the soft transition at the edge of the protection circle. Prevents hard "
        "edges between protected and reduced areas.\n\n"
        "Guide values: 6–10. For large, bright stars (long focal length, halos) rather 10–15."
    ),
    "search": (
        "Search radius (px)",
        "The catalogue says where a star should be. The script searches within this many "
        "pixels for the actual star. This compensates small errors from plate solving and "
        "lens distortion.\n\n"
        "Increase (30–40) if red circles appear at the image edge in the check image or the "
        "log reports a large position offset.\n"
        "Too large is bad in dense fields: a neighbouring star may be found instead.\n\n"
        "Focal length: wide angle with strong distortion rather 25–40, tele 10–20."
    ),
    "thresh": (
        "Edge threshold",
        "Brightness (0–1) at which a star ends for the radius measurement. Measured relative "
        "to the local background of the star image.\n\n"
        "Smaller (0.01–0.02): larger protection circles, more of the glow is kept.\n"
        "Larger (0.05–0.1): more compact protection circles.\n\n"
        "Works together with the halo transition: a pixel is only fully protected above "
        "edge threshold × halo transition."
    ),
    "halo": (
        "Halo transition",
        "Only the star itself is protected: the script measures the radial profile of each "
        "catalogue star. Neighbouring stars inside the protection area are reduced normally.\n\n"
        "The halo transition defines how much of the star's outer glow is kept:\n"
        "• profile above edge threshold × value: fully protected (star core)\n"
        "• below: smoothly reduced\n"
        "• below the edge threshold: reduced like the background\n\n"
        "With the defaults (0.03 × 6) everything above 0.18 is fully protected.\n\n"
        "Smaller (2–4): more glow around bright stars is kept.\n"
        "Larger (8–10): more compact stars, definitely no rings.\n"
        "1 = hard circle as in version 1.0."
    ),
    "glow": (
        "Glow strength",
        "Adds a soft glow around the stars – for the classic, slightly \"shining\" look of "
        "bright stars. 0 = off.\n\n"
        "For each catalogue star a round glow made of three Gaussian profiles is created – "
        "in the star's colour, sized to the star radius. Neighbouring stars do not create "
        "their own glow. It is blended with screen and soft saturation.\n\n"
        "Guide values: 0.5 subtle, 1–1.5 clear, 2+ strong.\n\n"
        "Two ways:\n"
        "• Apply: the glow is computed together with the reduction.\n"
        "• Apply glow only: adds the glow to the current image afterwards without reducing "
        "again. Can be repeated, Ctrl+Z undoes one step.\n\n"
        "Which stars get a glow is set by 'Glow from mag' and 'Glow to mag' – independent of "
        "the protection. If 'Glow catalogue stars only' is off, all bright stars in the image "
        "get a glow (the mag range is then ignored)."
    ),
    "glow_size": (
        "Glow size",
        "Size of the glow relative to the measured radius of each star. Brighter (larger) "
        "stars automatically get a larger glow – independent of focal length and resolution.\n\n"
        "Guide values: 1 = tight glow, 1.5–2 natural, 3+ distinct halo.\n\n"
        "With 'glow for all stars' the value is the base size in pixels (× 2)."
    ),
    "glow_soft": (
        "Glow softness",
        "Distributes the glow between core and falloff.\n\n"
        "• 0: compact, brighter glow close to the star\n"
        "• 0.5: balanced\n"
        "• 1: long, very soft falloff, like light haze\n\n"
        "For a natural look: softness 0.6–0.8 with moderate strength (0.5–1). The glow "
        "saturates softly and never clips hard at white."
    ),
    "glow_falloff": (
        "Glow falloff",
        "How quickly the glow gets weaker with decreasing brightness.\n\n"
        "Weight = ((Glow to − mag) / (Glow to − Glow from)) ^ falloff\n\n"
        "• at 'Glow from': full glow strength\n"
        "• at 'Glow to': zero\n\n"
        "Values:\n"
        "• 0: all stars in the range get the same glow\n"
        "• 0.5: drops slowly at first, then fast – medium stars shine too\n"
        "• 1: linear\n"
        "• 2–3: drops quickly – only the brightest stars shine clearly\n\n"
        "Example 'Glow from 0 to 6':\n"
        "a star of mag 3 gets half the strength at falloff 1, a quarter at 2, "
        "full strength at 0.\n\n"
        "Tip: set 'Glow from' to the brightness of the brightest star in the image "
        "(e.g. 0 instead of -2), otherwise no star reaches full strength.\n\n"
        "The size of the glow still follows the star radius."
    ),
    "glow_mag_min": (
        "Glow from mag",
        "Brightest stars that get a glow. Smaller number = brighter star.\n\n"
        "With -2 all bright stars are included. If you set e.g. 3, the very brightest stars "
        "(often already large and bright enough) get no additional glow.\n\n"
        "Examples:\n"
        "• protect to 6, glow -2 to 6: constellation stars shine\n"
        "• protect to 6, glow -2 to 8: medium stars get a delicate glow too\n"
        "• protect to 6, glow 6 to 9: only the fainter catalogue stars get a glow"
    ),
    "glow_mag_max": (
        "Glow to mag",
        "Faintest stars that get a glow – independent of the protection limit.\n\n"
        "How strongly the stars between 'Glow from' and 'Glow to' shine is controlled by "
        "the glow falloff.\n\n"
        "The catalogue is queried up to the larger of both limits (protection / glow). "
        "Hipparcos is only complete to about 7.5–8 mag – beyond that single stars are "
        "missing and the glow becomes patchy. Values above 8 also lengthen the star search."
    ),
    "preview": (
        "Preview",
        "Opens a preview window with an overview image and a crop as before/after.\n\n"
        "• Clicking into the overview chooses the crop (yellow frame).\n"
        "• Zoom: 1:1 shows real pixels, 2:1 enlarges for details, 1:2 shows more area.\n"
        "• Change parameters in the main window → the preview refreshes automatically "
        "after half a second.\n\n"
        "The star search runs once for the whole image and is cached. Only changes to "
        "magnitude, search radius or edge threshold trigger it again. Everything else is "
        "computed for the crop only and is therefore fast.\n\n"
        "The image in Siril is only changed by 'Apply'. 'Apply' always starts from the "
        "original – applying several times without Ctrl+Z is possible."
    ),
    "save": (
        "Save mask + check image",
        "Saves next to the image:\n\n"
        "• _protect_mask.fit – the effective protection mask (0–1), e.g. for your own "
        "PixelMath experiments\n"
        "• _protect_check.jpg – green = protected star found, red = not found\n"
        "• _protect_stars.csv – all catalogue stars with position, radius and status"
    ),
    "wcs": (
        "WCS from file",
        "The script needs the sky coordinates of the image (plate solve).\n\n"
        "Default: from the loaded image. If the star image is not solved, a plate-solved "
        "FITS of the same size can be chosen here, e.g. the original before StarNet."
    ),
    "focal": (
        "Presets",
        "Save parameter sets under any name, e.g. \"24mm Milky Way\", \"85mm Cygnus\" "
        "or \"Fiddlesticks\".\n\n"
        "• Drop-down: choose a saved preset – it is loaded immediately.\n"
        "• Save: save the current values under the typed name (existing ones are "
        "overwritten after confirmation).\n"
        "• Delete: remove the selected preset.\n"
        "• Factory defaults: the script's initial values.\n\n"
        "On start the last used preset is loaded. If there is none, a preset whose name "
        "contains the focal length of the image is searched (e.g. \"24mm\" or \"24 mm\").\n\n"
        "The focal length is read from the image (FOCALLEN) or estimated from pixel size and "
        "plate solve and is shown for information only.\n\n"
        "Stored in: %APPDATA%\\siril-starbeautifier\\presets.json"
    ),
    "lang": (
        "Language",
        "Switches the user interface between German and English. The choice is saved "
        "as default and used on the next start."
    ),
}

HELP_LANG_DE = (
    "Sprache",
    "Schaltet die Oberfläche zwischen Deutsch und Englisch um. Die Auswahl wird als "
    "Standard gespeichert und beim nächsten Start verwendet."
)
from astropy.coordinates import SkyCoord  # noqa: E402
import astropy.units as u  # noqa: E402

VERSION = "3.1"


# =============================================================================
# Kernlogik (ohne GUI, testbar)
# =============================================================================

def to_float_chw(data):
    """Pixeldaten -> float32 (C, H, W) im Bereich 0..1."""
    if np.issubdtype(data.dtype, np.integer):
        data = data.astype(np.float32) / 65535.0
    else:
        data = data.astype(np.float32)
        mx = float(np.nanmax(data)) if data.size else 1.0
        if mx > 1.0:
            data /= 65535.0 if mx <= 65535.0 else mx
    data = np.nan_to_num(np.clip(data, 0.0, 1.0))
    if data.ndim == 2:
        data = data[np.newaxis]
    return data


def wcs_from_header_string(hstr):
    if not hstr:
        return None
    sep = "\n" if "\n" in hstr else ""
    try:
        hdr = fits.Header.fromstring(hstr, sep=sep)
        w = WCS(hdr).celestial
    except Exception:
        return None
    if not w.has_celestial or not w.wcs.ctype[0]:
        return None
    return w


def wcs_from_file(path):
    with fits.open(path) as hdul:
        hdu = next((h for h in hdul if h.data is not None), hdul[0])
        hdr = hdu.header.copy()
        shape = hdu.data.shape[-2:] if hdu.data is not None else None
    w = WCS(hdr).celestial
    if not w.has_celestial or not w.wcs.ctype[0]:
        return None, shape
    return w, shape


def field_geometry(w, width, height):
    cx, cy = (width - 1) / 2.0, (height - 1) / 2.0
    center = w.pixel_to_world(cx, cy)
    xs = [0, width - 1, 0, width - 1, cx, cx, 0, width - 1]
    ys = [0, 0, height - 1, height - 1, 0, height - 1, cy, cy]
    edge = w.pixel_to_world(xs, ys)
    radius = float(center.separation(edge).deg.max()) + 1.0
    return center, radius


def query_hipparcos(center, radius_deg, maglim):
    from astroquery.vizier import Vizier
    v = Vizier(columns=["HIP", "_RAJ2000", "_DEJ2000", "Vmag"],
               column_filters={"Vmag": f"<{maglim}"}, row_limit=-1)
    res = v.query_region(center, radius=radius_deg * u.deg, catalog="I/239/hip_main")
    if len(res) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    t = res[0]
    vmag = np.asarray(t["Vmag"], dtype=float)
    ok = np.isfinite(vmag)
    return (np.asarray(t["_RAJ2000"], dtype=float)[ok],
            np.asarray(t["_DEJ2000"], dtype=float)[ok],
            vmag[ok], np.asarray(t["HIP"])[ok])


def load_catalog_csv(path):
    ra, dec, vmag, ids = [], [], [], []
    with open(path, newline="", encoding="utf-8") as fh:
        sample = fh.read(2048)
        fh.seek(0)
        delim = ";" if sample.count(";") >= sample.count(",") else ","
        for n, row in enumerate(csv.DictReader(fh, delimiter=delim)):
            row = {k.strip().lower(): v for k, v in row.items() if k}
            ra.append(float(row["ra"]))
            dec.append(float(row["dec"]))
            vmag.append(float(row["vmag"]))
            ids.append(int(row.get("id") or n))
    return np.array(ra), np.array(dec), np.array(vmag), np.array(ids)


def refine_and_measure(lum_s, lum, x0, y0, search, thresh, rmin, rmax):
    """
    Sucht das Maximum um (x0, y0), misst Radius und radiales Profil.
    -> (x, y, r, ok, profil) ; profil[k] = Median-Helligkeit im Ring k über dem
       lokalen Hintergrund (Nachbarsterne fallen durch den Median heraus)
    """
    h, w = lum.shape
    xi, yi = int(round(x0)), int(round(y0))
    x1, x2 = max(0, xi - search), min(w, xi + search + 1)
    y1, y2 = max(0, yi - search), min(h, yi + search + 1)
    if x1 >= x2 or y1 >= y2:
        return x0, y0, rmin, False, None
    box = lum_s[y1:y2, x1:x2]
    iy, ix = np.unravel_index(np.argmax(box), box.shape)
    px, py = x1 + ix, y1 + iy
    if box[iy, ix] < 2 * thresh:
        return x0, y0, rmin, False, None

    wx1, wx2 = max(0, px - rmax), min(w, px + rmax + 1)
    wy1, wy2 = max(0, py - rmax), min(h, py + rmax + 1)
    sub = lum[wy1:wy2, wx1:wx2]
    yy, xx = np.mgrid[wy1:wy2, wx1:wx2]
    d = np.hypot(xx - px, yy - py).astype(int)
    bg = float(np.median(sub))
    radius = rmax
    prof = [float(np.max(sub[d <= 1]))]
    for rr in range(1, rmax + 1):
        ring = sub[d == rr]
        if ring.size == 0:
            break
        med = float(np.median(ring))
        prof.append(med)
        if med < bg + thresh:
            radius = rr
            break
    # zwei Ringe Auslauf, damit das Profil weich auf 0 geht
    for rr in range(len(prof), len(prof) + 2):
        ring = sub[d == rr]
        prof.append(float(np.median(ring)) if ring.size else bg)
    prof = np.maximum(np.array(prof, np.float32) - bg, 0.0)
    prof[-1] = 0.0
    return float(px), float(py), float(max(rmin, radius)), True, prof


def locate_all(lum_s, lum, xs, ys, p):
    found = []
    for x0, y0 in zip(xs, ys):
        found.append(refine_and_measure(lum_s, lum, x0, y0, p["search"], p["thresh"],
                                        p["rmin"], p["rmax"]))
    return found


def match_quality(res, xs, ys, close=None):
    """
    -> (Anzahl, Median-Abweichung). Mit close werden nur Treffer gezählt, die
    näher als close px am Katalogort liegen (in dichten Feldern findet die
    Suche sonst in jeder Orientierung "irgendeinen" Stern).
    """
    ok = [(r[0] - x, r[1] - y) for r, x, y in zip(res, xs, ys) if r[3]]
    if not ok:
        return 0, 1e9
    sh = np.hypot(*np.array(ok).T)
    n = int((sh <= close).sum()) if close is not None else len(ok)
    return n, float(np.median(sh))


def draw_mask(shape, stars, margin, feather, ox=0, oy=0):
    """Schutzmaske für einen Bildbereich, der bei (ox, oy) beginnt."""
    h, w = shape
    mask = np.zeros(shape, np.float32)
    for x, y, r in stars:
        x, y = x - ox, y - oy
        R = r + margin + feather
        x1, x2 = int(max(0, x - R - 1)), int(min(w, x + R + 2))
        y1, y2 = int(max(0, y - R - 1)), int(min(h, y + R + 2))
        if x1 >= x2 or y1 >= y2:
            continue
        yy, xx = np.mgrid[y1:y2, x1:x2]
        d = np.hypot(xx - x, yy - y)
        m = np.clip((R - d) / max(feather, 1e-3), 0.0, 1.0)
        mask[y1:y2, x1:x2] = np.maximum(mask[y1:y2, x1:x2], m)
    return mask


def save_check(path, lum, found, missed, bottom_up, max_width=3000):
    from PIL import Image, ImageDraw
    disp = np.arcsinh(lum * 20) / np.arcsinh(20)
    img8 = (np.clip(disp, 0, 1) * 255).astype(np.uint8)
    h = img8.shape[0]
    if bottom_up:
        img8 = np.flipud(img8)
    im = Image.fromarray(img8).convert("RGB")
    sc = min(1.0, max_width / im.width)
    if sc < 1.0:
        im = im.resize((int(im.width * sc), int(im.height * sc)), Image.LANCZOS)
    dr = ImageDraw.Draw(im)
    for (x, y, r), col in [(f, (0, 255, 0)) for f in found] + [(m, (255, 0, 0)) for m in missed]:
        yy = (h - 1 - y) if bottom_up else y
        rr = max(3, r * sc)
        dr.ellipse([x * sc - rr, yy * sc - rr, x * sc + rr, yy * sc + rr], outline=col, width=2)
    im.save(path, quality=90)


def prepare(stars, w, p, log, catalog=None):
    """
    Katalog abfragen, Sterne im Bild finden und vermessen (ganzes Bild).
    Sucht alle Katalogsterne bis p["qmag"] (= max. aus Schutz- und Glow-Grenze).
    Hängt nur von qmag, search, thresh, rmin, rmax ab -> kann gecacht werden.
    -> dict mit found, missed, vmag_found, vmag_missed, profiles, rows, bottom_up, ...
    """
    qmag = p.get("qmag", p["maglim"])
    C, H, W = stars.shape
    center, radius = field_geometry(w, W, H)
    log(_("Feldmitte RA {ra:.2f}°, Dec {dec:.2f}°, Radius {rad:.1f}°",
          ra=center.ra.deg, dec=center.dec.deg, rad=radius))

    if catalog is None:
        log(_("Frage Hipparcos ab (V < {m}) ...", m=qmag))
        catalog = query_hipparcos(center, radius, qmag)
    ra, dec, vmag, ids = catalog
    sel = vmag < qmag
    ra, dec, vmag, ids = ra[sel], dec[sel], vmag[sel], ids[sel]

    coords = SkyCoord(ra * u.deg, dec * u.deg)
    in_cone = coords.separation(center).deg < radius
    px, py = w.world_to_pixel(coords)
    sr = p["search"]
    inside = in_cone & np.isfinite(px) & np.isfinite(py) & \
        (px > -sr) & (px < W + sr) & (py > -sr) & (py < H + sr)
    idx = np.nonzero(inside)[0]
    idx = idx[np.argsort(vmag[idx])]
    log(_("{n} Katalogsterne heller als {m} mag im Bild", n=len(idx), m=qmag))
    if qmag > 8.0:
        log(_("Hinweis: Hipparcos ist nur bis ca. 7.5–8 mag vollständig."))
    if len(idx) == 0:
        raise RuntimeError(_("Keine Katalogsterne im Bildfeld gefunden – Plate-Solve prüfen."))

    lum = stars.mean(axis=0)
    lum_s = gaussian_filter(lum, 1.5)

    # Zeilenorientierung automatisch bestimmen (FITS: Zeile 0 = unten),
    # nur mit den hellsten Sternen -> schnell auch bei Tausenden Katalogsternen
    xs = px[idx]
    no = min(len(idx), 300)
    ys_bu = py[idx]
    ys_td = (H - 1) - py[idx]
    close = max(3.0, 0.25 * p["search"])
    n_bu, med_bu = match_quality(locate_all(lum_s, lum, xs[:no], ys_bu[:no], p),
                                 xs[:no], ys_bu[:no], close)
    n_td, med_td = match_quality(locate_all(lum_s, lum, xs[:no], ys_td[:no], p),
                                 xs[:no], ys_td[:no], close)
    bottom_up = (n_bu, -med_bu) >= (n_td, -med_td)
    ys = ys_bu if bottom_up else ys_td
    log(_("Orientierung: {o} (enge Treffer {a} / {b} von {n}, Median {ma:.1f} / {mb:.1f} px)",
          o="bottom-up" if bottom_up else "top-down", a=n_bu, b=n_td, n=no, ma=med_bu, mb=med_td))
    res = locate_all(lum_s, lum, xs, ys, p)
    _n, med = match_quality(res, xs, ys)

    found, missed, rows = [], [], []
    prof_found, vm_found, vm_missed = [], [], []
    for k, i in enumerate(idx):
        x, y, r, ok, prof = res[k]
        if not ok:
            r = p["rmin"] + 2.0 * max(0.0, qmag - vmag[i])
        (found if ok else missed).append((x, y, r))
        (vm_found if ok else vm_missed).append(float(vmag[i]))
        if ok:
            prof_found.append(prof)
        rows.append([int(ids[i]), f"{vmag[i]:.2f}", f"{xs[k]:.1f}", f"{ys[k]:.1f}",
                     f"{x:.1f}", f"{y:.1f}", f"{r:.1f}", "ja" if ok else "nein"])
    log(_("Gefunden: {f}, nicht gefunden: {m}, Positionsabweichung Median {d:.1f} px",
          f=len(found), m=len(missed), d=med))
    if med > sr * 0.6:
        log(_("WARNUNG: große Positionsabweichung – Plate-Solve (Verzeichnung/SIP) prüfen "
              "oder Suchradius erhöhen."))

    return dict(found=found, missed=missed, rows=rows, bottom_up=bottom_up,
                catalog=catalog, lum=lum, profiles=prof_found + [None] * len(missed),
                vmag_found=np.array(vm_found), vmag_missed=np.array(vm_missed))


def protect_stars(prep, p):
    """Sterne + Profile, die geschützt werden (V < maglim)."""
    vm = np.concatenate([prep.get("vmag_found", np.zeros(len(prep["found"]))),
                         prep.get("vmag_missed", np.zeros(len(prep["missed"])))])
    allst = prep["found"] + prep["missed"]
    profs = prep.get("profiles") or [None] * len(allst)
    keep = [k for k in range(len(allst)) if vm[k] < p["maglim"]]
    return [allst[k] for k in keep], [profs[k] for k in keep]


def glow_weight(m, lo, hi, gamma):
    """
    Glow-Gewicht nach Helligkeit: 1 bei 'Glow von' (lo), 0 bei 'Glow bis' (hi).
    gamma = 0: konstant, 1: linear, >1: fällt schneller ab (nur Hellste leuchten).
    """
    if hi <= lo or gamma <= 0:
        return 1.0
    t = float(np.clip((hi - m) / (hi - lo), 0.0, 1.0))
    return t ** gamma


def glow_stars(prep, p):
    """
    Gefundene Sterne im Glow-Magnitudenbereich [glow_mag_min, glow_mag_max)
    -> Liste (x, y, r, gewicht)
    """
    vm = prep.get("vmag_found", np.zeros(len(prep["found"])))
    lo, hi = p.get("glow_mag_min", -5.0), p.get("glow_mag_max", p["maglim"])
    gamma = p.get("glow_falloff", 1.0)
    out = []
    for (x, y, r), m in zip(prep["found"], vm):
        if lo <= m < hi:
            wgt = glow_weight(m, lo, hi, gamma)
            if wgt > 1e-4:
                out.append((x, y, r, wgt))
    return out


def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def render(stars, prep, p, ox=0, oy=0):
    """
    Reduktion für einen Bildbereich (oder das ganze Bild) rechnen.
    stars: (C, h, w) Ausschnitt ab (ox, oy). -> (result, prot)

    Schutz nach Sternprofil: Für jeden geschützten Stern wird nur der Anteil
    eines Pixels geschützt, den das (radiale) Profil dieses Sterns erklärt.
    Nachbarsterne und Hintergrund innerhalb des Schutzbereichs werden ganz
    normal reduziert -> keine "Inseln" um geschützte Sterne.
    """
    C, H, W = stars.shape
    lum = stars.mean(axis=0)
    eps = 1e-6
    lsafe = np.maximum(lum, eps)
    F = np.power(lsafe, p["power"] - 1.0) * p["dim"]   # Reduktionsfaktor je Pixel

    halo = p.get("halo", 1.0)
    lo, hi = p["thresh"], p["thresh"] * max(halo, 1.0)
    margin, feather = p["margin"], max(p["feather"], 1e-3)

    P = np.zeros((H, W), np.float32)  # geschützte Helligkeit
    stars_all, profiles = protect_stars(prep, p)
    for (x, y, r), prof in zip(stars_all, profiles):
        xr, yr = x - ox, y - oy
        R = r + margin + feather
        x1, x2 = int(max(0, xr - R - 1)), int(min(W, xr + R + 2))
        y1, y2 = int(max(0, yr - R - 1)), int(min(H, yr + R + 2))
        if x1 >= x2 or y1 >= y2:
            continue
        yy, xx = np.ogrid[y1:y2, x1:x2]
        d = np.sqrt((xx - xr) ** 2 + (yy - yr) ** 2)
        wd = np.clip((R - d) / feather, 0.0, 1.0)
        Lw = lum[y1:y2, x1:x2]
        if prof is not None:
            S = np.interp(d, np.arange(len(prof)), prof, right=0.0)
            wl = _smoothstep((S - lo) / (hi - lo)) if halo > 1.0 else 1.0
            Pw = np.minimum(Lw, S) * wd * wl
        else:  # Stern nicht gefunden: Fallback nach Pixelhelligkeit
            wl = _smoothstep((Lw - lo) / (hi - lo)) if halo > 1.0 else 1.0
            Pw = Lw * wd * wl
        P[y1:y2, x1:x2] = np.maximum(P[y1:y2, x1:x2], Pw)

    T = P + (lum - P) * F            # Ziel-Helligkeit
    scale = T / lsafe
    result = np.clip(stars * scale[np.newaxis], 0.0, 1.0).astype(np.float32)
    prot = np.clip(P / lsafe, 0.0, 1.0).astype(np.float32)
    return result, prot


def compute(stars, w, p, log, catalog=None, prep=None):
    """Komplettlauf: prepare + render für das ganze Bild."""
    if prep is None:
        prep = prepare(stars, w, p, log, catalog)
    result, prot = render(stars, prep, p)
    out = dict(prep)
    out.update(result=result, mask=prot, prep=prep)
    return out


def prep_key(p, w):
    return (p.get("qmag", p["maglim"]), int(p["search"]), round(p["thresh"], 4),
            tuple(np.round(w.wcs.crval, 3)))


def _glow_weights(soft):
    """Gewichte der drei Gauss-Stufen; soft 0..1 verschiebt Energie in den Auslauf."""
    soft = float(np.clip(soft, 0.0, 1.0))
    w = np.array([0.6 - 0.4 * soft, 0.3, 0.1 + 0.4 * soft])
    return w / w.sum(), (1.0, 2.2, 4.5)


def star_colors(color_src, stars):
    """Mittlere Farbe je Stern (normiert auf hellsten Kanal)."""
    C, H, W = color_src.shape
    cols = []
    for st in stars:
        x, y, r = st[:3]
        rr = max(2, int(r * 0.7))
        xi, yi = int(round(x)), int(round(y))
        x1, x2 = max(0, xi - rr), min(W, xi + rr + 1)
        y1, y2 = max(0, yi - rr), min(H, yi + rr + 1)
        if x1 >= x2 or y1 >= y2:
            cols.append(np.ones(C, np.float32))
            continue
        m = color_src[:, y1:y2, x1:x2].reshape(C, -1).mean(axis=1)
        mx = float(m.max())
        cols.append((m / mx).astype(np.float32) if mx > 1e-4 else np.ones(C, np.float32))
    return cols


def synth_glow(shape, stars, colors, strength, size, soft, ox=0, oy=0):
    """
    Synthetischer, runder Glow je Stern (Summe aus drei Gauss-Profilen),
    Größe relativ zum gemessenen Sternradius, in Sternfarbe.
    shape: (C, h, w) des Zielbereichs, der bei (ox, oy) beginnt.
    -> Glow-Ebene 0..1 (weich gesättigt)
    """
    C, H, W = shape
    g = np.zeros(shape, np.float32)
    if not stars or strength <= 0:
        return g
    wts, scales = _glow_weights(soft)
    for st, col in zip(stars, colors):
        x, y, r = st[:3]
        wgt = st[3] if len(st) > 3 else 1.0
        sig = float(np.clip(size * r / 3.0, 1.0, 60.0))
        R = int(np.ceil(3.0 * scales[-1] * sig))
        xr, yr = x - ox, y - oy
        x1, x2 = int(max(0, np.floor(xr - R))), int(min(W, np.ceil(xr + R + 1)))
        y1, y2 = int(max(0, np.floor(yr - R))), int(min(H, np.ceil(yr + R + 1)))
        if x1 >= x2 or y1 >= y2:
            continue
        yy, xx = np.ogrid[y1:y2, x1:x2]
        d2 = ((xx - xr) ** 2 + (yy - yr) ** 2).astype(np.float32)
        prof = np.zeros(d2.shape, np.float32)
        edge = 0.0
        for wt, sc in zip(wts, scales):
            ss = sc * sig
            prof += wt * np.exp(-d2 / (2.0 * ss * ss))
            edge += wt * np.exp(-(R * R) / (2.0 * ss * ss))
        # Profil so verschieben, dass es am Radius R exakt 0 ist -> rund, ohne Kante
        # (sonst entstehen im AutoStretch sichtbare Quadrate am Fensterrand)
        prof = np.clip((prof - edge) / max(1.0 - edge, 1e-6), 0.0, None)
        # Stärke nach Helligkeit (Glow-Abnahme), Größe nach Sternradius
        amp = strength * 0.35 * wgt
        g[:, y1:y2, x1:x2] += amp * prof[np.newaxis] * col[:, np.newaxis, np.newaxis]
    return 1.0 - np.exp(-g)


def blur_glow(img, strength, size, soft):
    """Glow für alle Sterne über Weichzeichnung der hellen Kerne."""
    lum = img.mean(axis=0)
    src = img * (lum * lum)[np.newaxis]  # betont Sternkerne, unterdrückt Grundrauschen
    wts, scales = _glow_weights(soft)
    base = max(1.0, 2.0 * size)
    g = np.zeros_like(img)
    for wt, sc in zip(wts, scales):
        for c in range(img.shape[0]):
            g[c] += wt * gaussian_filter(src[c], base * sc)
    return 1.0 - np.exp(-g * strength * 3.0)


def blur_pad(size):
    return int(np.ceil(3.0 * 4.5 * max(1.0, 2.0 * size)))


def screen(img, glow):
    return np.clip(1.0 - (1.0 - img) * (1.0 - glow), 0.0, 1.0).astype(np.float32)


def add_glow(img, stars, color_src, p, only_protected=True, ox=0, oy=0):
    """
    img: (C,h,w) Zielbereich ab (ox, oy); stars: Liste (x, y, r) in Vollbild-Koordinaten;
    color_src: Vollbild (C,H,W) für die Sternfarben.
    """
    strength, size, soft = p.get("glow", 0.0), p.get("glow_size", 1.5), p.get("glow_soft", 0.5)
    if strength <= 0:
        return img
    if only_protected:
        if stars is None:
            raise RuntimeError(_("Für 'Glow nur für Katalogsterne' wird die Sternliste eines "
                                 "vorherigen Laufs benötigt – zuerst 'Anwenden' oder 'Vorschau'."))
        cols = star_colors(color_src, stars)
        glow = synth_glow(img.shape, stars, cols, strength, size, soft, ox, oy)
    else:
        glow = blur_glow(img, strength, size, soft)
    return screen(img, glow)


def save_side_files(base, out, log):
    fits.PrimaryHDU(out["mask"].astype(np.float32)).writeto(f"{base}_protect_mask.fit", overwrite=True)
    save_check(f"{base}_protect_check.jpg", out["lum"], out["found"], out["missed"], out["bottom_up"])
    with open(f"{base}_protect_stars.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter=";")
        wr.writerow(["HIP", "Vmag", "x_kat", "y_kat", "x_bild", "y_bild", "radius_px", "gefunden"])
        wr.writerows(out["rows"])
    log(_("Gespeichert: {b}_protect_mask.fit / _check.jpg / _stars.csv", b=base))


# =============================================================================
# Voreinstellungen je Brennweite
# =============================================================================

def _config_base():
    return os.environ.get("APPDATA") or os.path.join(os.path.expanduser("~"), ".config")


def preset_path():
    return os.path.join(_config_base(), "siril-starbeautifier", "presets.json")


def _old_preset_paths():
    # frühere Namen: "StarBeautifier" (v3.0) und "Star Protect Reduce" (bis v2.x)
    return [os.path.join(_config_base(), d, "presets.json")
            for d in ("siril-starbeautyfier", "siril-star-protect")]


def load_presets(path=None):
    """
    -> {"presets": {name: params}, "last": name | None}
    Ältere Dateien (Presets je Brennweite) werden automatisch übernommen.
    """
    if path is None:
        path = preset_path()
        if not os.path.exists(path):
            for old in _old_preset_paths():  # Übernahme aus früheren Versionen
                if os.path.exists(old):
                    path = old
                    break
    try:
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
    except (OSError, ValueError):
        d = {}
    if "presets" not in d:  # Migration v1.2–1.8
        presets = {}
        for k, v in (d.get("focal") or {}).items():
            presets[f"{k} mm"] = v
        if d.get("default"):
            presets["Standard"] = d["default"]
        d = {"presets": presets, "last": None}
    d.setdefault("presets", {})
    d.setdefault("last", None)
    return d


def save_presets(d, path=None):
    path = path or preset_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=2, ensure_ascii=False)
    os.replace(tmp, path)


def suggest_preset(d, focal=None):
    """Zuletzt benutztes Preset, sonst eines, dessen Name die Brennweite enthält."""
    names = list(d.get("presets", {}))
    if d.get("last") in names:
        return d["last"]
    if focal:
        import re
        f = int(round(focal))
        pat = re.compile(rf"(?<!\d){f}\s*mm", re.IGNORECASE)
        for n in names:
            if pat.search(n):
                return n
    return None


def focal_from_header(hstr, w=None):
    """Brennweite in mm aus FOCALLEN, sonst aus Pixelgröße + Plate-Solve-Skala."""
    if not hstr:
        return None
    sep = "\n" if "\n" in hstr else ""
    try:
        hdr = fits.Header.fromstring(hstr, sep=sep)
    except Exception:
        return None
    f = hdr.get("FOCALLEN")
    try:
        if f and float(f) > 0:
            return float(f)
    except (TypeError, ValueError):
        pass
    pix = hdr.get("XPIXSZ")
    if w is None:
        w = wcs_from_header_string(hstr)
    try:
        if pix and w is not None:
            from astropy.wcs.utils import proj_plane_pixel_scales
            scale = float(np.mean(proj_plane_pixel_scales(w))) * 3600.0  # "/px
            if scale > 0:
                return 206.265 * float(pix) / scale
    except Exception:
        pass
    return None


# =============================================================================
# GUI
# =============================================================================

PARAMS = [
    # key, label, default, from, to, increment, tooltip
    ("maglim", "Schutz bis Magnitude", 6.0, 1.0, 9.0, 0.5,
     "Sterne heller als dieser V-Wert bleiben unverändert (Hipparcos, vollständig bis ca. 7.5)."),
    ("power", "Exponent (Stärke)", 2.0, 1.0, 5.0, 0.25,
     "Dämpfung der übrigen Sterne: L^Exponent. 2 = moderat, 3 = stark."),
    ("dim", "Helligkeitsfaktor", 1.0, 0.1, 1.0, 0.05,
     "Zusätzlicher Faktor für die übrigen Sterne."),
    ("margin", "Maskenrand (px)", 3.0, 0.0, 30.0, 1.0,
     "Vergrößert die Schutzkreise, falls Halos angeschnitten werden."),
    ("feather", "Weiche Kante (px)", 8.0, 0.0, 30.0, 1.0,
     "Breite des weichen Übergangs am Maskenrand."),
    ("search", "Suchradius (px)", 20, 3, 100, 1,
     "Wie weit vom Katalogort nach dem Stern gesucht wird (Restverzeichnung)."),
    ("thresh", "Randschwelle", 0.03, 0.005, 0.3, 0.005,
     "Helligkeit, ab der der Sternrand endet (für die Radiusmessung)."),
    ("halo", "Halo-Übergang", 6.0, 1.0, 20.0, 0.5,
     "Schutz folgt dem Sternprofil: unterhalb Randschwelle x Wert wird auch im Schutzkreis "
     "weich reduziert. Höher = weicherer Übergang, 1 = harter Kreis (alte Version)."),
]

GLOW_PARAMS = [
    ("glow", "Glow-Stärke", 0.0, 0.0, 3.0, 0.1,
     "Zusätzlicher weicher Schein um die Sterne. 0 = aus."),
    ("glow_size", "Glow-Größe", 1.5, 0.3, 6.0, 0.1,
     "Größe des Scheins relativ zum Sternradius."),
    ("glow_soft", "Glow-Weichheit", 0.5, 0.0, 1.0, 0.05,
     "0 = kompakter Schein, 1 = langer, weicher Auslauf."),
    ("glow_mag_min", "Glow von Mag", -2.0, -2.0, 12.0, 0.5,
     "Hellste Sterne, die einen Glow bekommen (kleinere Zahl = heller)."),
    ("glow_mag_max", "Glow bis Mag", 6.0, -1.0, 12.0, 0.5,
     "Schwächste Sterne, die einen Glow bekommen."),
    ("glow_falloff", "Glow-Abnahme", 1.0, 0.0, 5.0, 0.1,
     "0 = alle gleich, 1 = linear von 'Glow von' (voll) bis 'Glow bis' (null), >1 = schneller."),
]
ALL_PARAMS = PARAMS + GLOW_PARAMS

HELP_OVERVIEW = (
    "StarBeautifier\n\n"
    "Helle Sterne aus dem Hipparcos-Katalog werden geschützt, alle anderen Sterne "
    "im Sternbild werden gedämpft. So bleiben die Sternbilder erkennbar, während "
    "der dichte Sternteppich zurücktritt und Nebel und Dunkelwolken mehr Raum bekommen.\n\n"
    "Ablauf:\n"
    "1. Sternbild aus StarNet laden (nicht das Gesamtbild).\n"
    "2. Plate Solving mit Verzeichnungskorrektur (SIP).\n"
    "3. Werte einstellen, Anwenden, Kontrollbild prüfen.\n"
    "4. Mit Strg+Z zurück, Werte ändern, erneut anwenden.\n"
    "5. Danach in PixelMath: starless + stars\n\n"
    "Fahre mit der Maus über einen Parameter, um hier die Erklärung zu sehen."
)

HELP = {
    "maglim": (
        "Schutz bis Magnitude",
        "Sterne, die heller sind als dieser Wert (V-Magnitude laut Hipparcos), bleiben "
        "unverändert. Kleinere Zahl = hellerer Stern.\n\n"
        "Richtwerte:\n"
        "• 4–5: nur die markanten Sternbildsterne\n"
        "• 5.5–6: alles, was unter dunklem Himmel mit bloßem Auge sichtbar ist\n"
        "• 7–8: sehr viele geschützte Sterne – die Reduktion wirkt dann kaum noch "
        "und das Feld kann ungleichmäßig aussehen\n\n"
        "Brennweite: Bei Weitwinkel (14–24 mm) eher 5–5.5, weil sonst sehr viele "
        "Sterne im Feld geschützt werden. Bei Normal- und Teleobjektiven (50 mm+) "
        "ist das Feld kleiner, dort passen 6–7.\n\n"
        "Hipparcos ist bis ca. 7.5 mag vollständig."
    ),
    "power": (
        "Exponent (Stärke)",
        "Wie stark die nicht geschützten Sterne gedämpft werden. Gerechnet wird "
        "Helligkeit hoch Exponent – schwache Sterne verlieren dabei überproportional, "
        "helle nur wenig. Die Sternfarbe bleibt erhalten.\n\n"
        "Beispiel: Pixelhelligkeit 0.3\n"
        "• Exponent 1.5 → 0.16\n"
        "• Exponent 2 → 0.09\n"
        "• Exponent 3 → 0.03\n\n"
        "Richtwerte: 1.5 dezent, 2 moderat, 2.5–3 kräftig. Über 3 wirkt das Feld "
        "schnell leer und unnatürlich."
    ),
    "dim": (
        "Helligkeitsfaktor",
        "Zusätzlicher linearer Faktor für alle nicht geschützten Sterne, nach dem "
        "Exponenten. 1.0 = keine zusätzliche Abdunklung.\n\n"
        "Unterschied zum Exponenten: Der Faktor dämpft alle Sterne gleichmäßig, der "
        "Exponent vor allem die schwachen. Mit 0.7–0.8 lassen sich auch die mittelhellen "
        "Sterne zurücknehmen, ohne den Exponenten weiter zu erhöhen."
    ),
    "margin": (
        "Maskenrand (px)",
        "Vergrößert den gemessenen Schutzkreis jedes hellen Sterns um diese Pixelzahl.\n\n"
        "Erhöhen, wenn die Ränder heller Sterne angeschnitten oder ausgefranst wirken. "
        "Seit dem Halo-Übergang ist ein großer Rand unkritisch, weil schwache Pixel im "
        "Kreis trotzdem reduziert werden.\n\n"
        "Brennweite: Mit längerer Brennweite und höherer Auflösung werden Sterne in "
        "Pixeln größer – dann eher 5–10 statt 3."
    ),
    "feather": (
        "Weiche Kante (px)",
        "Breite des weichen Übergangs am Rand des Schutzkreises. Verhindert harte "
        "Kanten zwischen geschütztem und reduziertem Bereich.\n\n"
        "Richtwerte: 6–10. Bei großen, hellen Sternen (lange Brennweite, Halos) "
        "eher 10–15."
    ),
    "search": (
        "Suchradius (px)",
        "Der Katalog sagt, wo ein Stern sein müsste. Das Script sucht im Umkreis von "
        "so vielen Pixeln nach dem tatsächlichen Stern. Das gleicht kleine Fehler aus "
        "Plate-Solving und Objektivverzeichnung aus.\n\n"
        "Erhöhen (30–40), wenn im Kontrollbild rote Kreise am Bildrand auftauchen "
        "oder das Protokoll eine große Positionsabweichung meldet.\n"
        "Zu groß ist schlecht in dichten Feldern: Dann wird evtl. ein Nachbarstern "
        "gefunden.\n\n"
        "Brennweite: Weitwinkel mit starker Verzeichnung eher 25–40, Tele 10–20."
    ),
    "thresh": (
        "Randschwelle",
        "Helligkeit (0–1), bei der ein Stern für die Radiusmessung endet. Gemessen "
        "wird relativ zum lokalen Hintergrund des Sternbilds.\n\n"
        "Kleiner (0.01–0.02): größere Schutzkreise, mehr vom Schein bleibt erhalten.\n"
        "Größer (0.05–0.1): kompaktere Schutzkreise.\n\n"
        "Wirkt zusammen mit dem Halo-Übergang: Erst ab Randschwelle × Halo-Übergang "
        "ist ein Pixel voll geschützt."
    ),
    "halo": (
        "Halo-Übergang",
        "Geschützt wird nur der Stern selbst: Das Script misst das radiale Profil "
        "jedes Katalogsterns. Nachbarsterne im Schutzbereich werden normal reduziert.\n\n"
        "Der Halo-Übergang bestimmt, wie viel vom äußeren Schein des Sterns erhalten bleibt:\n"
        "• Profil über Randschwelle × Wert: voll geschützt (Sternkern)\n"
        "• darunter: fließend reduziert\n"
        "• unter der Randschwelle: wie der Hintergrund reduziert\n\n"
        "Mit den Standardwerten (0.03 × 6) ist alles über 0.18 voll geschützt.\n\n"
        "Kleiner (2–4): mehr Schein um helle Sterne bleibt erhalten.\n"
        "Größer (8–10): kompaktere Sterne, sicher keine Ringe.\n"
        "1 = harter Kreis wie in Version 1.0."
    ),
    "glow": (
        "Glow-Stärke",
        "Legt einen weichen Schein um die Sterne – für den klassischen, leicht "
        "\"leuchtenden\" Look heller Sterne. 0 = aus.\n\n"
        "Für jeden geschützten Stern wird ein runder Schein aus drei Gauss-Profilen "
        "erzeugt – in der Farbe des Sterns, mit Größe passend zum Sternradius. "
        "Nachbarsterne erzeugen keinen eigenen Schein. Überlagert wird per Screen "
        "mit sanfter Sättigung.\n\n"
        "Richtwerte: 0.5 dezent, 1–1.5 deutlich, 2+ kräftig.\n\n"
        "Zwei Wege:\n"
        "• Anwenden: Glow wird direkt mit der Reduktion gerechnet.\n"
        "• Nur Glow anwenden: legt den Glow nachträglich auf das aktuelle Bild, ohne "
        "erneut zu reduzieren. Mehrfach möglich, Strg+Z macht einen Schritt rückgängig.\n\n"
        "Welche Sterne einen Glow bekommen, legen 'Glow von Mag' und 'Glow bis Mag' "
        "fest – unabhängig vom Schutz. Ist 'Glow nur für Katalogsterne' aus, bekommen "
        "alle hellen Sterne im Bild einen Schein (Mag-Bereich wird dann ignoriert)."
    ),
    "glow_size": (
        "Glow-Größe",
        "Größe des Scheins relativ zum gemessenen Radius jedes Sterns. Hellere (größere) "
        "Sterne bekommen dadurch automatisch einen größeren Schein – unabhängig von "
        "Brennweite und Auflösung.\n\n"
        "Richtwerte: 1 = enger Schein, 1.5–2 natürlich, 3+ deutlicher Hof.\n\n"
        "Bei 'Glow für alle Sterne' gilt der Wert als Grundgröße in Pixeln (× 2)."
    ),
    "glow_soft": (
        "Glow-Weichheit",
        "Verteilt den Schein zwischen Kern und Auslauf.\n\n"
        "• 0: kompakter, hellerer Schein dicht am Stern\n"
        "• 0.5: ausgewogen\n"
        "• 1: langer, sehr weicher Auslauf, wie leichter Dunst\n\n"
        "Für einen natürlichen Look: Weichheit 0.6–0.8 mit moderater Stärke (0.5–1). "
        "Der Schein sättigt sanft und schneidet nie hart bei Weiß ab."
    ),
    "glow_falloff": (
        "Glow-Abnahme",
        "Wie stark der Glow mit abnehmender Helligkeit schwächer wird.\n\n"
        "Gewicht = ((Glow bis − Mag) / (Glow bis − Glow von)) ^ Abnahme\n\n"
        "• bei 'Glow von': volle Glow-Stärke\n"
        "• bei 'Glow bis': null\n\n"
        "Werte:\n"
        "• 0: alle Sterne im Bereich bekommen gleich viel Glow\n"
        "• 0.5: fällt erst langsam, dann schnell ab – auch mittlere Sterne leuchten\n"
        "• 1: linear\n"
        "• 2–3: fällt schnell ab – nur die hellsten Sterne leuchten deutlich\n\n"
        "Beispiel 'Glow von 0, bis 6':\n"
        "Stern mit Mag 3 bekommt bei Abnahme 1 die halbe, bei 2 ein Viertel, "
        "bei 0 die volle Stärke.\n\n"
        "Tipp: 'Glow von' auf die Helligkeit des hellsten Sterns im Bild setzen "
        "(z. B. 0 statt -2), sonst erreicht kein Stern die volle Stärke.\n\n"
        "Die Größe des Scheins richtet sich weiterhin nach dem Sternradius."
    ),
    "glow_mag_min": (
        "Glow von Mag",
        "Hellste Sterne, die einen Glow bekommen. Kleinere Zahl = hellerer Stern.\n\n"
        "Mit -2 sind alle hellen Sterne dabei. Setzt du z. B. 3, bekommen die "
        "allerhellsten Sterne (die oft schon groß und hell genug sind) keinen "
        "zusätzlichen Schein.\n\n"
        "Beispiele:\n"
        "• Schutz bis 6, Glow -2 bis 6: Sternbildsterne leuchten\n"
        "• Schutz bis 6, Glow -2 bis 8: auch mittlere Sterne bekommen einen zarten Schein\n"
        "• Schutz bis 6, Glow 6 bis 9: nur die schwächeren Katalogsterne bekommen Glow"
    ),
    "glow_mag_max": (
        "Glow bis Mag",
        "Schwächste Sterne, die einen Glow bekommen – unabhängig von der Schutzgrenze.\n\n"
        "Wie stark die Sterne zwischen 'Glow von' und 'Glow bis' leuchten, steuert "
        "die Glow-Abnahme.\n\n"
        "Der Katalog wird bis zur größeren der beiden Grenzen (Schutz / Glow) "
        "abgefragt. Hipparcos ist nur bis ca. 7.5–8 mag vollständig – darüber "
        "fehlen einzelne Sterne, der Glow wird dann lückenhaft. Werte über 8 "
        "verlängern außerdem die Sternsuche."
    ),
    "preview": (
        "Vorschau",
        "Öffnet ein Vorschaufenster mit Übersichtsbild und einem Ausschnitt in "
        "Vorher/Nachher.\n\n"
        "• Klick ins Übersichtsbild wählt den Ausschnitt (gelber Rahmen).\n"
        "• Zoom: 1:1 zeigt echte Pixel, 2:1 vergrößert für Details, 1:2 zeigt mehr Fläche.\n"
        "• Parameter im Hauptfenster ändern → die Vorschau aktualisiert sich nach "
        "einer halben Sekunde automatisch.\n\n"
        "Die Sternsuche läuft einmal für das ganze Bild und wird zwischengespeichert. "
        "Nur Änderungen an Magnitude, Suchradius oder Randschwelle lösen sie erneut aus. "
        "Alles andere wird nur für den Ausschnitt gerechnet und ist deshalb schnell.\n\n"
        "Das Bild in Siril wird erst mit 'Anwenden' verändert. 'Anwenden' rechnet "
        "immer vom Original – mehrfaches Anwenden ohne Strg+Z ist möglich."
    ),
    "save": (
        "Maske + Kontrollbild speichern",
        "Speichert neben dem Bild:\n\n"
        "• _protect_mask.fit – die wirksame Schutzmaske (0–1), z. B. für eigene "
        "PixelMath-Versuche\n"
        "• _protect_check.jpg – grün = geschützter Stern gefunden, rot = nicht gefunden\n"
        "• _protect_stars.csv – alle Katalogsterne mit Position, Radius und Status"
    ),
    "wcs": (
        "WCS aus Datei",
        "Das Script braucht die Himmelskoordinaten des Bildes (Plate Solve).\n\n"
        "Standard: aus dem geladenen Bild. Ist das Sternbild nicht gesolvt, kann hier "
        "ein geplatesolvtes FITS gleicher Größe gewählt werden, z. B. das Original vor "
        "StarNet."
    ),
    "focal": (
        "Voreinstellungen",
        "Parametersätze unter frei wählbarem Namen speichern, z. B. "
        "\"24mm Milchstraße\", \"85mm Cygnus\" oder \"Pustekuchen\".\n\n"
        "• Auswahlliste: gespeichertes Preset wählen – es wird sofort geladen.\n"
        "• Speichern: aktuelle Werte unter dem eingetippten Namen speichern "
        "(vorhandene werden nach Rückfrage überschrieben).\n"
        "• Löschen: gewähltes Preset entfernen.\n"
        "• Werkseinstellungen: Ausgangswerte des Scripts.\n\n"
        "Beim Start wird das zuletzt benutzte Preset geladen. Gibt es keins, wird ein "
        "Preset gesucht, dessen Name die Brennweite des Bildes enthält (z. B. \"24mm\" "
        "oder \"24 mm\").\n\n"
        "Die Brennweite wird aus dem Bild gelesen (FOCALLEN) oder aus Pixelgröße und "
        "Plate-Solve geschätzt und nur zur Information angezeigt.\n\n"
        "Gespeichert in: %APPDATA%\\siril-starbeautifier\\presets.json"
    ),
}


HELP["lang"] = HELP_LANG_DE


def help_text(key):
    if key is None:
        return (_("Übersicht"), HELP_OVERVIEW_EN if LANG == "en" else HELP_OVERVIEW)
    return (HELP_EN if LANG == "en" else HELP)[key]


class ParamControl:
    """Schieberegler + Eingabefeld (statt Spinbox), Mausrad = ein Schritt."""

    def __init__(self, parent, default, lo, hi, inc):
        self.lo, self.hi, self.inc = float(lo), float(hi), float(inc)
        istr = repr(float(inc))
        self.dec = 0 if float(inc).is_integer() else len(istr.split(".")[1])
        self.var = tk.DoubleVar(value=self._snap(default))
        self.frame = ttk.Frame(parent)
        self.scale = ttk.Scale(self.frame, from_=self.lo, to=self.hi, orient="horizontal",
                               length=130, command=self._on_scale)
        self.scale.pack(side=tk.LEFT, padx=(0, 6))
        self.text = tk.StringVar(value=self._fmt(self.var.get()))
        self.entry = ttk.Entry(self.frame, textvariable=self.text, width=7, justify="right")
        self.entry.pack(side=tk.LEFT)
        self._busy = False
        self.scale.set(self.var.get())
        self.entry.bind("<Return>", self._on_entry)
        self.entry.bind("<FocusOut>", self._on_entry)
        for w in (self.entry, self.scale):
            w.bind("<MouseWheel>", self._on_wheel)                       # Windows / macOS
            w.bind("<Button-4>", lambda e: self._step(+1))               # Linux
            w.bind("<Button-5>", lambda e: self._step(-1))
        self.entry.bind("<Up>", lambda e: self._step(+1))
        self.entry.bind("<Down>", lambda e: self._step(-1))
        self.var.trace_add("write", self._on_var)

    def _snap(self, v):
        v = min(max(float(v), self.lo), self.hi)
        v = self.lo + round((v - self.lo) / self.inc) * self.inc
        return round(min(max(v, self.lo), self.hi), max(self.dec, 0) + 2)

    def _fmt(self, v):
        txt = f"{v:.4f}".rstrip("0").rstrip(".")
        if "." not in txt and self.dec > 0:
            txt += "." + "0" * self.dec
        elif "." in txt and len(txt.split(".")[1]) < self.dec:
            txt += "0" * (self.dec - len(txt.split(".")[1]))
        return txt

    def _on_scale(self, val):
        if self._busy:
            return
        v = self._snap(val)
        if abs(v - self.var.get()) > 1e-9:
            self.var.set(v)

    def _on_entry(self, _e=None):
        try:
            v = float(self.text.get().replace(",", "."))
            v = round(min(max(v, self.lo), self.hi), 4)  # getippte Werte nur begrenzen
        except ValueError:
            v = self.var.get()
        self.var.set(v)
        self.text.set(self._fmt(v))

    def _step(self, n):
        self.var.set(self._snap(self.var.get() + n * self.inc))
        return "break"

    def _on_wheel(self, e):
        return self._step(+1 if e.delta > 0 else -1)

    def _on_var(self, *_):
        try:
            v = float(self.var.get())
        except (tk.TclError, ValueError):
            return
        self._busy = True
        try:
            self.scale.set(v)
        finally:
            self._busy = False
        self.text.set(self._fmt(v))


class StarProtectGUI:
    def __init__(self, root):
        self.root = root
        self.root.title(f"StarBeautifier v{VERSION}")
        self.root.resizable(False, False)
        self.wcs_file = None
        self.style = tksiril.standard_style()

        self.ok = False
        self.siril = s.SirilInterface()
        try:
            self.siril.connect()
        except s.SirilConnectionError:
            messagebox.showerror(_("Fehler"), _("Keine Verbindung zu Siril."))
            self.root.destroy()
            return
        try:
            self.siril.cmd("requires", "1.4.0")
        except s.CommandError:
            self.root.destroy()
            return
        if not self.siril.is_image_loaded():
            self.siril.error_messagebox(_("Kein Bild geladen. Bitte zuerst das Sternbild öffnen."))
            self.root.destroy()
            return

        self.wcs_file = None
        self.catalog_cache = {}  # (maglim, center) -> catalog
        self.vars = {}
        self.source = None
        self.last_result = None
        self.last_prot = None
        self.prep_cache = {}
        self.preview = None
        self._preview_job = None
        self.presets = load_presets()
        global LANG
        LANG = self.presets.get("language") or detect_language()
        self._build()
        self._init_focal_and_preset()
        self._add_traces()
        tksiril.match_theme_to_siril(self.root, self.siril)
        self.ok = True

    # ---------------------------------------------------------------- UI
    def _build(self):
        outer = ttk.Frame(self.root, padding=10)
        outer.pack(fill=tk.BOTH, expand=True)
        frm = ttk.Frame(outer)
        frm.grid(row=0, column=0, sticky="nsew")
        self._build_help_panel(outer)

        ttk.Label(frm, text=_("Helle Katalogsterne schützen, Rest reduzieren")).grid(
            row=0, column=0, sticky="w", pady=(0, 8))
        lf = ttk.Frame(frm)
        lf.grid(row=0, column=1, sticky="e", pady=(0, 8))
        ll = ttk.Label(lf, text=_("Sprache") + ":")
        ll.pack(side=tk.LEFT, padx=(0, 4))
        self.lang_var = tk.StringVar(value="English" if LANG == "en" else "Deutsch")
        lcb = ttk.Combobox(lf, textvariable=self.lang_var, values=["Deutsch", "English"],
                           state="readonly", width=9)
        lcb.pack(side=tk.LEFT)
        lcb.bind("<<ComboboxSelected>>", lambda _e: self._switch_language())
        self._bind_help("lang", ll, lcb)

        def add_rows(params, r):
            for key, label, default, lo, hi, inc, tip in params:
                lab = ttk.Label(frm, text=_(label))
                lab.grid(row=r, column=0, sticky="w", pady=2)
                ctl = ParamControl(frm, default, lo, hi, inc)
                ctl.frame.grid(row=r, column=1, sticky="e", pady=2)
                tksiril.create_tooltip(ctl.entry, _(tip))
                self.vars[key] = ctl.var
                self._bind_help(key, lab, ctl.scale, ctl.entry)
                r += 1
            return r

        r = add_rows(PARAMS, 1)
        ttk.Separator(frm).grid(row=r, column=0, columnspan=2, sticky="ew", pady=(8, 4))
        r += 1
        gl = ttk.Label(frm, text="Glow")
        gl.grid(row=r, column=0, sticky="w")
        self._bind_help("glow", gl)
        r = add_rows(GLOW_PARAMS, r + 1)
        self.glow_prot_var = tk.BooleanVar(value=True)
        gcb = ttk.Checkbutton(frm, text=_("Glow nur für Katalogsterne (Mag-Bereich)"), variable=self.glow_prot_var)
        gcb.grid(row=r, column=0, sticky="w", pady=2)
        self.glow_btn = ttk.Button(frm, text=_("Nur Glow anwenden"), command=self._apply_glow_only)
        self.glow_btn.grid(row=r, column=1, sticky="e", pady=2)
        self._bind_help("glow", gcb, self.glow_btn)

        r += 1
        self.save_var = tk.BooleanVar(value=True)
        cb = ttk.Checkbutton(frm, text=_("Maske + Kontrollbild speichern"), variable=self.save_var)
        cb.grid(row=r, column=0, columnspan=2, sticky="w", pady=(8, 2))
        self._bind_help("save", cb)

        r += 1
        self.wcs_label = ttk.Label(frm, text=_("WCS: aus geladenem Bild") if not self.wcs_file
                                   else f"WCS: {os.path.basename(self.wcs_file)}")
        self.wcs_label.grid(row=r, column=0, sticky="w", pady=2)
        b = ttk.Button(frm, text=_("WCS aus Datei…"), command=self._pick_wcs)
        b.grid(row=r, column=1, sticky="e", pady=2)
        self._bind_help("wcs", b, self.wcs_label)
        tksiril.create_tooltip(b, _("Nur nötig, wenn das geladene Bild nicht plate-gesolvt ist: "
                                    "geplatesolvtes FITS gleicher Größe wählen."))

        r += 1
        ttk.Separator(frm).grid(row=r, column=0, columnspan=2, sticky="ew", pady=8)
        r += 1
        pl = ttk.Label(frm, text=_("Voreinstellung"))
        pl.grid(row=r, column=0, sticky="w", pady=2)
        self.preset_name = tk.StringVar(value="")
        self.preset_cb = ttk.Combobox(frm, textvariable=self.preset_name, width=22)
        self.preset_cb.grid(row=r, column=1, sticky="e", pady=2)
        self.preset_cb.bind("<<ComboboxSelected>>", lambda _e: self._load_preset())
        tksiril.create_tooltip(self.preset_cb, _("Preset wählen oder neuen Namen eintippen."))
        r += 1
        pf = ttk.Frame(frm)
        pf.grid(row=r, column=0, columnspan=2, sticky="ew", pady=2)
        b1 = ttk.Button(pf, text=_("Speichern"), command=self._save_preset)
        b1.pack(side=tk.LEFT, padx=(0, 4))
        b2 = ttk.Button(pf, text=_("Löschen"), command=self._delete_preset)
        b2.pack(side=tk.LEFT, padx=4)
        b3 = ttk.Button(pf, text=_("Werkseinstellungen"), command=self._factory_reset)
        b3.pack(side=tk.LEFT, padx=4)
        r += 1
        self.focal_var = tk.DoubleVar(value=0.0)
        self.preset_label = ttk.Label(frm, text=_("Brennweite laut Bild: unbekannt"))
        self.preset_label.grid(row=r, column=0, columnspan=2, sticky="w", pady=2)
        self._bind_help("focal", pl, self.preset_cb, b1, b2, b3, self.preset_label)

        r += 1
        bf = ttk.Frame(frm)
        bf.grid(row=r, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(bf, text=_("Schließen"), command=self._close).pack(side=tk.RIGHT, padx=4)
        pv = ttk.Button(bf, text=_("Vorschau…"), command=self._open_preview)
        pv.pack(side=tk.LEFT, padx=4)
        self._bind_help("preview", pv)
        self.apply_btn = ttk.Button(bf, text=_("Anwenden"), command=self._apply)
        self.apply_btn.pack(side=tk.RIGHT, padx=4)

    # ------------------------------------------------------- Bild & Cache
    def _params(self):
        p = {k: float(v.get()) for k, v in self.vars.items()}
        p["search"] = int(round(p["search"]))
        p["rmin"], p["rmax"] = 4, 80
        p["qmag"] = max(p["maglim"], p.get("glow_mag_max", p["maglim"]))
        return p

    def _get_source(self):
        """
        Liefert das Original-Sternbild. Ist das aktuelle Bild noch unser letztes
        Ergebnis, wird das gespeicherte Original verwendet (kein Strg+Z nötig).
        """
        data = np.asarray(self.siril.get_image_pixeldata())
        cur = to_float_chw(data)
        if self.source is not None and self.last_result is not None \
                and cur.shape == self.last_result.shape:
            a = cur[:, ::37, ::37]
            b = self.last_result[:, ::37, ::37]
            if np.allclose(a, b, atol=2e-3):
                return self.source, data.ndim
        if self.source is None or cur.shape != self.source.shape \
                or not np.allclose(cur[:, ::37, ::37], self.source[:, ::37, ::37], atol=2e-3):
            self.prep_cache = {}
            self.last_prot = None
        self.source = cur
        self.last_result = None
        return cur, data.ndim

    def _get_wcs(self, H, W):
        if self.wcs_file:
            w, shape = wcs_from_file(self.wcs_file)
            if w is None:
                raise RuntimeError(_("Die gewählte Datei enthält kein WCS."))
            if shape is not None and tuple(shape) != (H, W):
                raise RuntimeError(_("Größe passt nicht: WCS-Datei {a}x{b}, Bild {c}x{d}",
                                     a=shape[1], b=shape[0], c=W, d=H))
        else:
            hstr = self.siril.get_image_fits_header(return_as="str")
            cached = getattr(self, "_wcs_cache", None)
            if cached is not None and cached[0] == hstr:
                w = cached[1]
            else:
                w = wcs_from_header_string(hstr)
                self._wcs_cache = (hstr, w)
            if w is None:
                raise RuntimeError(_("Das Bild ist nicht plate-gesolvt. Bitte zuerst Plate Solving "
                                     "ausführen oder 'WCS aus Datei…' wählen."))
        if not self._focal():
            try:
                f = focal_from_header(self.siril.get_image_fits_header(return_as="str"), w)
            except Exception:
                f = None
            if f:
                self._set_focal(f)
                self._log(_("StarBeautifier: Brennweite aus Plate-Solve geschätzt: {f:.1f} mm", f=f))
        return w

    def _get_prep(self, stars, w, p):
        key = prep_key(p, w)
        prep = self.prep_cache.get(key)
        if prep is None:
            self._log(_("StarBeautifier: WCS gelesen, SIP-Verzeichnung: {s}",
                        s=_("ja") if w.sip is not None else _("nein")))
            ckey = (p["qmag"], tuple(np.round(w.wcs.crval, 3)))
            prep = prepare(stars, w, p, self._log, catalog=self.catalog_cache.get(ckey))
            self.catalog_cache[ckey] = prep["catalog"]
            self.prep_cache[key] = prep
        self.last_prep = prep
        return prep

    def _set_image(self, result, ndim, undo_label):
        self.last_result = result
        out = result[0] if ndim == 2 else result
        with self.siril.image_lock():
            try:
                self.siril.undo_save_state(undo_label)
            except Exception:
                pass
            self.siril.set_image_pixeldata(out)

    # ------------------------------------------------------------ Vorschau
    def _open_preview(self):
        if self.preview is not None and self.preview.alive():
            self.preview.refresh(full=True)
            self.preview.lift()
            return
        try:
            stars, _nd = self._get_source()
            C, H, W = stars.shape
            w = self._get_wcs(H, W)
            p = self._params()
            prep = self._get_prep(stars, w, p)
        except Exception as e:
            self._log(_("StarBeautifier Fehler: {e}", e=e))
            messagebox.showerror("StarBeautifier", str(e))
            return
        self.preview = PreviewWindow(self, stars, prep["bottom_up"])
        self.preview.refresh(full=True)

    def _preview_render(self, x1, y1, x2, y2):
        """Vorher/Nachher für den Ausschnitt [y1:y2, x1:x2] (Array-Koordinaten)."""
        stars = self.source
        C, H, W = stars.shape
        p = self._params()
        w = self._get_wcs(H, W)
        prep = self._get_prep(stars, w, p)
        glow_on = p.get("glow", 0) > 0
        only_prot = bool(self.glow_prot_var.get())
        pad = blur_pad(p["glow_size"]) if (glow_on and not only_prot) else 0
        px1, py1 = max(0, x1 - pad), max(0, y1 - pad)
        px2, py2 = min(W, x2 + pad), min(H, y2 + pad)
        crop = stars[:, py1:py2, px1:px2]
        res, prot = render(crop, prep, p, px1, py1)
        if glow_on:
            res = add_glow(res, glow_stars(prep, p), stars, p, only_prot, px1, py1)
        ys, xs = slice(y1 - py1, y2 - py1), slice(x1 - px1, x2 - px1)
        return stars[:, y1:y2, x1:x2], res[:, ys, xs], prep["bottom_up"]

    def _on_param_change(self, *_):
        if self.preview is not None and self.preview.alive():
            if self._preview_job is not None:
                self.root.after_cancel(self._preview_job)
            self._preview_job = self.root.after(500, self._debounced_preview)

    def _debounced_preview(self):
        self._preview_job = None
        if self.preview is not None and self.preview.alive():
            self.preview.refresh()

    def _apply_glow_only(self):
        self.glow_btn.state(["disabled"])
        self.root.update()
        try:
            p = self._params()
            strength = p["glow"]
            if strength <= 0:
                raise RuntimeError(_("Glow-Stärke ist 0 – bitte einen Wert > 0 einstellen."))
            data = self.siril.get_image_pixeldata()
            img = to_float_chw(np.asarray(data))
            only_prot = bool(self.glow_prot_var.get())
            prep = getattr(self, "last_prep", None)
            if only_prot and (self.source is None or self.source.shape != img.shape):
                prep = None
            self._progress(_("StarBeautifier: Glow …"), 0.4)
            res = add_glow(img, glow_stars(prep, p) if prep else None,
                           self.source if prep else img, p, only_prot)
            self._set_image(res, np.asarray(data).ndim, "StarBeautifier Glow")
            self._log(_("StarBeautifier: Glow {g} (Größe {s}, Weichheit {w}) nachträglich angewendet",
                        g=strength, s=p["glow_size"], w=p["glow_soft"]))
        except Exception as e:
            self._log(_("StarBeautifier Fehler: {e}", e=e))
            messagebox.showerror("StarBeautifier", str(e))
        finally:
            try:
                self.siril.reset_progress()
            except Exception:
                pass
            try:
                self.glow_btn.state(["!disabled"])
            except Exception:
                pass

    # --------------------------------------------------------- Hilfe-Panel
    def _build_help_panel(self, parent):
        box = ttk.LabelFrame(parent, text=_("Erklärung"), padding=6)
        box.grid(row=0, column=1, sticky="nsew", padx=(12, 0))
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        style = ttk.Style()
        bg = style.lookup("TFrame", "background") or "#2b2b2b"
        fg = style.lookup("TLabel", "foreground") or "#e0e0e0"
        import tkinter.font as tkfont
        try:
            base_font = tkfont.nametofont("TkDefaultFont")
            self._title_font = base_font.copy()
            self._title_font.configure(weight="bold", size=abs(int(base_font.cget("size"))) + 1)
        except Exception:
            base_font, self._title_font = None, None
        self.help_text = tk.Text(box, wrap="word", width=46, height=30, relief="flat",
                                 font=base_font, spacing1=1, spacing3=1,
                                 borderwidth=0, padx=4, pady=4, bg=bg, fg=fg,
                                 cursor="arrow", highlightthickness=0)
        sb = ttk.Scrollbar(box, orient="vertical", command=self.help_text.yview)
        self.help_text.configure(yscrollcommand=sb.set)
        self.help_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        try:
            if self._title_font is not None:
                self.help_text.tag_configure("title", font=self._title_font, spacing3=8)
        except Exception:
            pass
        self._show_help(None)

    def _show_help(self, key):
        title, body = help_text(key)
        t = self.help_text
        t.configure(state="normal")
        t.delete("1.0", "end")
        t.insert("end", title + "\n", "title")
        t.insert("end", body)
        t.configure(state="disabled")

    def _bind_help(self, key, *widgets):
        for wdg in widgets:
            for ev in ("<Enter>", "<FocusIn>"):
                wdg.bind(ev, lambda _e, k=key: self._show_help(k), add="+")

    def _add_traces(self):
        for v in list(self.vars.values()) + [self.glow_prot_var]:
            v.trace_add("write", self._on_param_change)

    def _switch_language(self):
        global LANG
        new = "en" if self.lang_var.get() == "English" else "de"
        if new == LANG:
            return
        params = self._current_params()
        name = self.preset_name.get()
        focal = self._focal()
        if self.preview is not None and self.preview.alive():
            self.preview.close()
        self.preview = None
        LANG = new
        # als Standard speichern
        try:
            d = load_presets()
            d["language"] = LANG
            save_presets(d)
            self.presets = d
        except OSError:
            pass
        for child in self.root.winfo_children():
            child.destroy()
        self.vars = {}
        self._build()
        self._set_params(params)
        self._refresh_preset_list()
        self.preset_name.set(name)
        if focal:
            self._set_focal(focal)
        self._add_traces()
        try:
            tksiril.match_theme_to_siril(self.root, self.siril)
        except Exception:
            pass

    # ------------------------------------------------------------ Presets
    def _current_params(self):
        d = {k: round(float(v.get()), 4) for k, v in self.vars.items()}
        d["save_files"] = bool(self.save_var.get())
        d["glow_only_protected"] = bool(self.glow_prot_var.get())
        return d

    def _set_params(self, d):
        for key, _l, default, *_ in ALL_PARAMS:
            self.vars[key].set(d.get(key, default) if d else default)
        if d and "save_files" in d:
            self.save_var.set(bool(d["save_files"]))
        if d and "glow_only_protected" in d:
            self.glow_prot_var.set(bool(d["glow_only_protected"]))

    def _focal(self):
        try:
            f = float(self.focal_var.get())
        except Exception:
            f = 0.0
        return f if f > 0 else None

    def _set_focal(self, f):
        self.focal_var.set(round(f, 1))
        try:
            self.preset_label.config(text=_("Brennweite laut Bild: {f:.0f} mm", f=f))
        except Exception:
            pass

    def _refresh_preset_list(self):
        names = sorted(self.presets["presets"], key=str.lower)
        try:
            self.preset_cb["values"] = names
        except Exception:
            pass

    def _init_focal_and_preset(self):
        f = None
        try:
            f = focal_from_header(self.siril.get_image_fits_header(return_as="str"))
        except Exception:
            pass
        if f:
            self._set_focal(f)
            self._log(_("StarBeautifier: Brennweite aus Bild: {f:.1f} mm", f=f))
        self.presets = load_presets()
        self._refresh_preset_list()
        name = suggest_preset(self.presets, f)
        if name:
            self.preset_name.set(name)
            self._load_preset(quiet=True)

    def _load_preset(self, quiet=False):
        self.presets = load_presets()
        name = self.preset_name.get().strip()
        d = self.presets["presets"].get(name)
        if d is None:
            if not quiet:
                messagebox.showwarning("StarBeautifier", _("Kein Preset namens „{n}“.", n=name))
            return
        self._set_params(d)
        self.presets["last"] = name
        try:
            save_presets(self.presets)
        except OSError:
            pass
        self._log(_("StarBeautifier: Preset „{n}“ geladen", n=name))

    def _save_preset(self):
        name = self.preset_name.get().strip()
        if not name:
            messagebox.showwarning("StarBeautifier",
                                   _("Bitte einen Namen für das Preset eintippen."))
            return
        self.presets = load_presets()
        if name in self.presets["presets"] and not messagebox.askyesno(
                "StarBeautifier", _("Preset „{n}“ überschreiben?", n=name)):
            return
        self.presets["presets"][name] = self._current_params()
        self.presets["last"] = name
        save_presets(self.presets)
        self._refresh_preset_list()
        self._log(_("StarBeautifier: Preset „{n}“ gespeichert ({p})", n=name, p=preset_path()))

    def _delete_preset(self):
        name = self.preset_name.get().strip()
        self.presets = load_presets()
        if name not in self.presets["presets"]:
            messagebox.showwarning("StarBeautifier", _("Kein Preset namens „{n}“.", n=name))
            return
        if not messagebox.askyesno("StarBeautifier", _("Preset „{n}“ löschen?", n=name)):
            return
        del self.presets["presets"][name]
        if self.presets.get("last") == name:
            self.presets["last"] = None
        save_presets(self.presets)
        self._refresh_preset_list()
        self.preset_name.set("")
        self._log(_("StarBeautifier: Preset „{n}“ gelöscht", n=name))

    def _factory_reset(self):
        self._set_params(None)
        self.save_var.set(True)
        self.glow_prot_var.set(True)
        self._log(_("StarBeautifier: Werkseinstellungen geladen"))

    def _pick_wcs(self):
        path = filedialog.askopenfilename(title="Geplatesolvtes FITS wählen",
                                          filetypes=[("FITS", "*.fit *.fits *.fts"), ("Alle", "*.*")])
        if path:
            self.wcs_file = path
            self.wcs_label.config(text=f"WCS: {os.path.basename(path)}")

    def _close(self):
        try:
            self.siril.disconnect()
        except Exception:
            pass
        self.root.quit()
        self.root.destroy()

    def _log(self, msg):
        try:
            self.siril.log(msg)
        except Exception:
            print(msg)

    def _progress(self, msg, frac):
        try:
            self.siril.update_progress(msg, frac)
        except Exception:
            pass
        self.root.update_idletasks()

    # ------------------------------------------------------------- Apply
    def _apply(self):
        self.apply_btn.state(["disabled"])
        self.root.update()
        try:
            p = self._params()
            self._progress(_("StarBeautifier: lese Bild …"), 0.05)
            stars, ndim = self._get_source()
            C, H, W = stars.shape
            w = self._get_wcs(H, W)
            self._progress(_("StarBeautifier: Katalog & Sternsuche …"), 0.3)
            prep = self._get_prep(stars, w, p)
            self._progress(_("StarBeautifier: Reduktion …"), 0.6)
            out = compute(stars, w, p, self._log, prep=prep)

            self.last_prot = out["mask"]
            self.last_prep = prep
            self.last_shape = (H, W)
            result = out["result"]
            if p.get("glow", 0) > 0:
                self._progress(_("StarBeautifier: Glow …"), 0.75)
                result = add_glow(result, glow_stars(prep, p), stars, p,
                                  bool(self.glow_prot_var.get()))
                self._log(_("StarBeautifier: Glow {g} (Größe {s}, Weichheit {w}) angewendet",
                            g=p["glow"], s=p["glow_size"], w=p["glow_soft"]))
            self._progress(_("StarBeautifier: schreibe Ergebnis …"), 0.85)
            self._set_image(result, ndim, "StarBeautifier")

            if self.save_var.get():
                fname = None
                try:
                    fname = self.siril.get_image_filename()
                except Exception:
                    pass
                if fname:
                    base = os.path.splitext(fname)[0]
                else:
                    base = os.path.join(os.getcwd(), "star_protect")
                save_side_files(base, out, self._log)

            self._log(_("StarBeautifier: fertig. Erneutes 'Anwenden' rechnet wieder vom Original."))
        except (SirilError, RuntimeError, Exception) as e:  # noqa: B014
            self._log(_("StarBeautifier Fehler: {e}", e=e))
            messagebox.showerror("StarBeautifier", str(e))
        finally:
            try:
                self.siril.reset_progress()
            except Exception:
                pass
            self.apply_btn.state(["!disabled"])


class PreviewWindow:
    PANEL = 480
    ZOOMS = (("1:1", 1.0), ("2:1 (Zoom)", 2.0), ("1:2 (Übersicht)", 0.5))

    def _zoom_map(self):
        return {_(k): v for k, v in self.ZOOMS}

    def __init__(self, gui, stars, bottom_up):
        from PIL import Image, ImageTk  # noqa: F401
        self.gui = gui
        self.H, self.W = stars.shape[1:]
        self.bottom_up = bottom_up
        self.center = (self.W // 2, self.H // 2)  # Array-Koordinaten
        self.top = tk.Toplevel(gui.root)
        self.top.title(_("StarBeautifier – Vorschau"))
        self.top.protocol("WM_DELETE_WINDOW", self.close)
        self._alive = True
        self._imgs = {}

        frm = ttk.Frame(self.top, padding=8)
        frm.pack(fill=tk.BOTH, expand=True)

        # Übersicht
        ov = self._to_display(stars)
        ov_im = Image.fromarray(ov)
        sc = min(self.PANEL * 2 / ov_im.width, 260 / ov_im.height, 1.0)
        self.ov_scale = sc
        ov_im = ov_im.resize((max(1, int(ov_im.width * sc)), max(1, int(ov_im.height * sc))),
                             Image.LANCZOS)
        self._imgs["ov"] = ImageTk.PhotoImage(ov_im)
        top_row = ttk.Frame(frm)
        top_row.grid(row=0, column=0, columnspan=2, sticky="w")
        self.ov_canvas = tk.Canvas(top_row, width=ov_im.width, height=ov_im.height,
                                   highlightthickness=0, cursor="crosshair")
        self.ov_canvas.create_image(0, 0, anchor="nw", image=self._imgs["ov"])
        self.ov_canvas.pack(side=tk.LEFT)
        self.ov_canvas.bind("<Button-1>", self._on_click)
        self.rect = None

        ctl = ttk.Frame(top_row, padding=(10, 0))
        ctl.pack(side=tk.LEFT, anchor="n")
        ttk.Label(ctl, text=_("Klick ins Bild wählt\nden Ausschnitt.")).pack(anchor="w", pady=(0, 8))
        ttk.Label(ctl, text="Zoom").pack(anchor="w")
        self.zoom_var = tk.StringVar(value=_("1:1"))
        cb = ttk.Combobox(ctl, textvariable=self.zoom_var, values=list(self._zoom_map()),
                          state="readonly", width=14)
        cb.pack(anchor="w", pady=(0, 8))
        cb.bind("<<ComboboxSelected>>", lambda _e: self.refresh(full=True))
        ttk.Button(ctl, text=_("Aktualisieren"), command=lambda: self.refresh(full=True))\
            .pack(anchor="w", pady=2)
        self.status = ttk.Label(ctl, text="")
        self.status.pack(anchor="w", pady=(8, 0))

        # Vorher / Nachher
        ttk.Label(frm, text=_("Vorher")).grid(row=1, column=0, sticky="w", pady=(8, 2))
        ttk.Label(frm, text=_("Nachher")).grid(row=1, column=1, sticky="w", pady=(8, 2), padx=(8, 0))
        self.before = tk.Canvas(frm, width=self.PANEL, height=self.PANEL, highlightthickness=0, bg="black")
        self.after = tk.Canvas(frm, width=self.PANEL, height=self.PANEL, highlightthickness=0, bg="black")
        self.before.grid(row=2, column=0)
        self.after.grid(row=2, column=1, padx=(8, 0))
        try:
            tksiril.match_theme_to_siril(self.top, gui.siril)
        except Exception:
            pass

    def alive(self):
        return self._alive

    def lift(self):
        self.top.lift()

    def close(self):
        self._alive = False
        self.top.destroy()

    def _to_display(self, arr):
        """(C,h,w) 0..1 -> uint8 (h,w,3) in Siril-Ansicht (oben = oben)."""
        a = np.clip(arr, 0, 1)
        if a.shape[0] == 1:
            a = np.repeat(a, 3, axis=0)
        img = (np.transpose(a[:3], (1, 2, 0)) * 255).astype(np.uint8)
        return np.flipud(img) if self.bottom_up else img

    def _on_click(self, ev):
        x = int(ev.x / self.ov_scale)
        row_disp = int(ev.y / self.ov_scale)
        y = (self.H - 1 - row_disp) if self.bottom_up else row_disp
        self.center = (min(max(x, 0), self.W - 1), min(max(y, 0), self.H - 1))
        self.refresh()

    def _crop_box(self):
        z = self._zoom_map().get(self.zoom_var.get(), 1.0)
        size = int(round(self.PANEL / z))
        size = min(size, self.W, self.H)
        cx, cy = self.center
        x1 = int(min(max(cx - size // 2, 0), self.W - size))
        y1 = int(min(max(cy - size // 2, 0), self.H - size))
        return x1, y1, x1 + size, y1 + size

    def _draw_rect(self, x1, y1, x2, y2):
        sc = self.ov_scale
        if self.bottom_up:
            dy1, dy2 = self.H - y2, self.H - y1
        else:
            dy1, dy2 = y1, y2
        if self.rect is not None:
            self.ov_canvas.delete(self.rect)
        self.rect = self.ov_canvas.create_rectangle(x1 * sc, dy1 * sc, x2 * sc, dy2 * sc,
                                                    outline="#ffd400", width=2)

    def refresh(self, full=False):
        from PIL import Image, ImageTk
        import time
        if not self._alive:
            return
        x1, y1, x2, y2 = self._crop_box()
        self._draw_rect(x1, y1, x2, y2)
        self.status.config(text=_("rechne …"))
        self.top.update_idletasks()
        t0 = time.time()
        try:
            before, after, self.bottom_up = self.gui._preview_render(x1, y1, x2, y2)
        except Exception as e:
            self.status.config(text=_("Fehler"))
            self.gui._log(_("StarBeautifier Vorschau-Fehler: {e}", e=e))
            return
        for key, arr, canvas in (("b", before, self.before), ("a", after, self.after)):
            im = Image.fromarray(self._to_display(arr))
            if im.width != self.PANEL:
                rs = Image.NEAREST if im.width < self.PANEL else Image.LANCZOS
                im = im.resize((self.PANEL, self.PANEL), rs)
            self._imgs[key] = ImageTk.PhotoImage(im)
            canvas.delete("all")
            canvas.create_image(0, 0, anchor="nw", image=self._imgs[key])
        self.status.config(text=f"{(time.time() - t0) * 1000:.0f} ms")


def main():
    try:
        root = ThemedTk()
        gui = StarProtectGUI(root)
        if gui.ok:
            root.mainloop()
    except SirilError as e:
        print(f"Fehler: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
