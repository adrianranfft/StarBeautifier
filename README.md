# StarBeautifier

*by Adrian Ranfft*

A Siril script that tames dense star fields in wide-field and Milky Way images while keeping the constellations intact – and optionally adds a soft, natural glow to the bright stars.

Bright stars are identified by **real astronomical magnitude** using the Hipparcos catalogue and the plate solve of your image. Everything fainter is reduced, so nebulae and dark clouds get more room while Vega, Deneb, Antares & co. stay exactly as they are.

---

## Features

- **Catalogue-based star protection** – stars brighter than a chosen magnitude (Hipparcos, via VizieR) are protected; all others are reduced.
- **Profile-based protection** – only the protected star itself is kept. Faint neighbouring stars right next to it are reduced normally, so there are no "islands" or rings around protected stars.
- **Colour-preserving reduction** – stars are dimmed via their luminance (`L^exponent`), star colours stay intact.
- **Synthetic star glow** – a round, soft glow per star in the star's own colour, sized to the measured star radius, with its own magnitude range and a brightness-dependent falloff.
- **Live preview** – before/after view of a selectable crop, updating automatically while you move the sliders.
- **Presets** – save parameter sets under any name (e.g. "16mm Milky Way", "85mm Cygnus").
- **German / English user interface**, switchable at any time.
- **Non-destructive workflow** – "Apply" always starts from the original star image; Ctrl+Z works as usual.
- Check image, protection mask and star list (CSV) are saved for inspection.

---

## Requirements

- **Siril 1.4 or newer** with Python script support (included in the standard Windows, macOS and Linux builds).
- **Internet connection**
  - on first start, to install the required Python packages,
  - for the Hipparcos catalogue query (VizieR, once per magnitude limit and session).
- **A separated star image**, e.g. from StarNet in Siril, StarXTerminator or similar. The script works on the *stars only* image, not on the combined image.
- **A plate-solved image.** The star image itself can be plate-solved in Siril (a stars-only image solves well). Enable distortion correction (SIP) for wide-angle lenses.

The following Python packages are installed automatically by Siril on first start:
`numpy`, `scipy`, `astropy`, `astroquery`, `pillow`, `ttkthemes`.

---

## Installation

1. Download `StarBeautifier.py` from this repository.
2. Put it into a folder of your choice, e.g. `Documents/Siril-Scripts`.
3. In Siril open **Preferences → Scripts** and add this folder to the script directories.
4. Restart Siril. The script now appears in the **Scripts** menu.

On the first start, Siril installs the required Python packages. This can take a minute – watch the Siril log.

**Updating:** simply replace the file. A restart of Siril is not necessary; just close the script window and start it again from the menu.

---

## Usage

### Typical workflow

1. **Separate stars and background**
   Run StarNet (or another star removal tool) on your stretched image. You get a *starless* image and a *star* image (`starmask_…`).
2. **Open the star image** in Siril.
3. **Plate solve it**: *Tools → Astrometry → Image Plate Solver*, with distortion correction (SIP) enabled.
   Alternatively, use the button **WCS from file…** in the script to take the plate solve from another FITS file of the same size (e.g. the original image before StarNet).
4. **Start the script**: *Scripts → StarBeautifier*.
5. **Adjust the parameters**. Use **Preview…** to see the effect on a crop before applying.
6. Click **Apply**. The star image in Siril is replaced by the result.
   You can change parameters and click **Apply** again – it always starts from the original star image.
7. **Recombine** with the starless image, e.g. in PixelMath:
   ```
   starless + stars
   ```
   or with Siril's star recomposition tool.

### Preview window

- Click into the overview image to choose the crop (yellow frame).
- Zoom: **1:1** real pixels, **2:1** enlarged, **1:2** larger area.
- Changing any parameter in the main window updates the preview after half a second.
- The full star search runs only once and is cached; only changes to magnitude limits, search radius or edge threshold trigger a new search.

### Recommended display mode

Judge the result with Siril's display set to **Linear**. After a strong reduction the star image's background is almost exactly zero, and **AutoStretch** then amplifies even tiny values to white, which can make the image look very different from the actual data.

---

## Parameters

Hover over any parameter in the script window – the explanation panel on the right shows a detailed description.

### Star protection & reduction

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Protect up to magnitude** | 6.0 | 1 – 9 | Stars brighter than this V magnitude (Hipparcos) stay unchanged. Smaller number = brighter star. 4–4.5 keeps the constellation figures, 5.5–6 everything visible to the naked eye. Values of 7–8 protect so many stars that the reduction barely works. |
| **Exponent (strength)** | 2.0 | 1 – 5 | Reduction of all unprotected stars: `L^exponent`. Faint stars lose disproportionately, bright ones only a little. 1.5 subtle, 2 moderate, 2.5–3 strong. Example: pixel brightness 0.3 → 0.09 at exponent 2. |
| **Brightness factor** | 1.0 | 0.1 – 1 | Additional linear dimming of all unprotected stars, applied after the exponent. Useful to tone down medium-bright stars without raising the exponent. |
| **Mask margin (px)** | 3 | 0 – 30 | Enlarges the protection area around each protected star. Increase if star edges look clipped. |
| **Feather (px)** | 8 | 0 – 30 | Width of the soft transition at the edge of the protection area. |
| **Search radius (px)** | 20 | 3 – 100 | How far from the catalogue position the script searches for the real star (compensates plate solve and distortion residuals). Wide angle: 25–40, tele: 10–20. Too large values can pick up neighbouring stars in dense fields. |
| **Edge threshold** | 0.03 | 0.005 – 0.3 | Brightness (relative to the local background) at which a star ends for the radius measurement. Smaller = larger protected area, more glow kept. |
| **Halo transition** | 6.0 | 1 – 20 | How much of a protected star's outer glow is kept. The star's profile is fully protected above *edge threshold × value* and smoothly reduced below. Smaller (2–4) keeps more glow, larger (8–10) gives compact stars. 1 = hard circle. |

### Glow

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Glow strength** | 0.0 | 0 – 3 | Strength of the added glow. 0 = off. 0.5 subtle, 1–1.5 clear, 2+ strong. |
| **Glow size** | 1.5 | 0.3 – 6 | Size of the glow relative to the measured star radius – brighter (larger) stars automatically get a larger glow. 1 tight, 1.5–2 natural, 3+ distinct halo. |
| **Glow softness** | 0.5 | 0 – 1 | Distribution between core and falloff. 0 = compact, bright glow; 1 = long, very soft falloff like light haze. |
| **Glow from mag** | -2.0 | -2 – 12 | Brightest stars that get a glow; they receive the full glow strength. Tip: set this to the brightest star in your image (e.g. 0), otherwise no star reaches full strength. |
| **Glow to mag** | 6.0 | -1 – 12 | Faintest stars that get a glow; the glow reaches zero here. Independent of the protection limit. |
| **Glow falloff** | 1.0 | 0 – 5 | How the glow strength decreases from *Glow from* (full) to *Glow to* (zero): `weight = ((to − mag) / (to − from)) ^ falloff`. 0 = all stars equal, 1 = linear, 2–3 = only the brightest stars shine clearly. |

**Glow catalogue stars only (mag range)** – when enabled (default), only catalogue stars within the glow magnitude range get a glow. When disabled, all bright stars in the image get a glow based on a blur of their cores, and the magnitude range is ignored.

**Apply glow only** – adds the glow to the current image without running the reduction again. Can be repeated; Ctrl+Z undoes one step.

### Other options

| Option | Description |
|---|---|
| **Save mask + check image** | Saves the files listed below next to the image. |
| **WCS from file…** | Use the plate solve of another FITS file of the same size, if the loaded star image is not plate-solved. |
| **Preset** | Choose a saved preset (loads immediately) or type a new name. **Save**, **Delete**, **Factory defaults**. On start, the last used preset is loaded; if there is none, a preset whose name contains the image's focal length (e.g. "16mm") is suggested. |
| **Language** | Switches between German and English. The choice is saved as default. |

### Starting points by focal length

| Lens | Protect up to | Exponent | Search radius | Glow |
|---|---|---|---|---|
| 14–24 mm (Milky Way, wide field) | 4.5 – 5.5 | 2 – 2.5 | 25 – 40 | from 0 to 4, falloff 1–2 |
| 35–50 mm | 5.5 – 6 | 2 | 20 – 25 | from 0 to 5, falloff 1 |
| 85 mm and longer | 6 – 7 | 1.75 – 2 | 10 – 20 | from 2 to 6, falloff 1 |

These are only starting points – use the preview and save what works as a preset.

---

## Output files

When *Save mask + check image* is enabled, the following files are written next to the image:

| File | Content |
|---|---|
| `<name>_protect_mask.fit` | Effective protection mask (0–1), e.g. for your own PixelMath experiments |
| `<name>_protect_check.jpg` | Check image: green circle = protected star found, red = not found |
| `<name>_protect_stars.csv` | All catalogue stars with HIP number, magnitude, catalogue and measured position, radius and status |

Presets and the language setting are stored in:

- Windows: `%APPDATA%\siril-starbeautifier\presets.json`
- macOS / Linux: `~/.config/siril-starbeautifier/presets.json`

---

## How it works

1. The WCS (plate solve) is read from the loaded image.
2. Hipparcos stars inside the field are queried from VizieR (catalogue `I/239/hip_main`).
3. Each catalogue star is located in the star image (local search around the predicted position) and its radius and radial profile are measured. Row orientation (bottom-up / top-down) is detected automatically.
4. Reduction: every pixel's luminance is reduced to `L^exponent × factor`, except the part of a pixel that is explained by the radial profile of a protected star.
5. Glow: for each star in the glow range, three Gaussian profiles in the star's colour are combined, weighted by magnitude, softly saturated and blended with *screen*.

---

## Troubleshooting

**"The image is not plate-solved"**
Plate solve the star image in Siril first, or choose a plate-solved FITS of the same size via *WCS from file…*.

**Red circles in the check image**
Those catalogue stars were not found near their predicted position. Stars below the horizon or outside the sharp image area are normal. Red circles at the image edges usually mean distortion residuals: increase the *search radius* or enable SIP distortion correction for plate solving. A "large position offset" warning in the log points to the same problem.

**Many protected stars still look patchy or unnatural**
The protection limit is probably too faint (7–8 mag). Use 4.5–6.

**Glow looks foggy**
Too many stars with too much glow overlap. Lower *Glow to mag*, raise *Glow falloff*, or reduce *Glow size*.

**Strange squares or bright areas after Apply**
Switch Siril's display mode from *AutoStretch* to *Linear* – see [Recommended display mode](#recommended-display-mode).

**Catalogue query fails**
The Hipparcos data are fetched online from VizieR (CDS Strasbourg). Check your internet connection or firewall and try again.

**Note on faint limits**
Hipparcos is complete to about 7.5–8 mag. Limits beyond 8 are possible but the star list becomes incomplete and the star search takes longer.

---

## Author

**Adrian Ranfft**

---

## Credits

- [Siril](https://siril.org) and its Python interface `sirilpy`
- Hipparcos catalogue (ESA 1997), accessed via [VizieR](https://vizier.cds.unistra.fr), CDS Strasbourg
- [Astropy](https://www.astropy.org), [astroquery](https://astroquery.readthedocs.io), NumPy, SciPy, Pillow, ttkthemes

---

## License

Copyright © 2026 Adrian Ranfft

This project is licensed under the **GNU General Public License v3.0** – see [LICENSE](LICENSE) for details.
