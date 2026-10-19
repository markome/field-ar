# Field Lines AR: browser-based Hiro-marker prototype

A small A-Frame + AR.js application that **overlays a precomputed magnetic-field-line GLB** on a marker. No native mobile app or electromagnetic calculation in the browser is required.

This is a **starter implementation**, not yet tested on your particular phone or laboratory apparatus. The bundled geometry is a *demo* of field lines around an ideal circular current loop, computed using the axisymmetric Biot–Savart field; it is not calibrated to your physical coil.

## Files

- `index.html` — browser-based AR scene with marker tracking and show/hide control
- `models/demo_fieldlines.glb` — illustrative field lines (tube meshes) and an orange reference loop
- `print-marker.html` — printable link to the official Hiro marker (loads image from AR.js)
- `tools/create_demo_model.py` — optional source for generating the demo GLB; uses NumPy, SciPy, Trimesh

## Quick start on a mobile phone

1. Extract the files and upload the **contents** of this directory into a GitHub repository or another static HTTPS host.
2. On GitHub: **Settings → Pages → Deploy from a branch → main → /(root) → Save**. Wait for the HTTPS URL.
3. Print or display the Hiro marker: <https://raw.githubusercontent.com/AR-js-org/AR.js/master/data/images/hiro.png>. The provided `print-marker.html` also shows it.
4. On a phone, open the HTTPS URL. Grant permission to use the camera, and point it at the marker.
5. A circular current loop and colored calculated magnetic field lines should appear above the marker. Press **Hide field** to compare with the real view.

**Important:** most mobile browsers only allow camera capture on HTTPS websites (or `localhost` on that very device). Visiting `http://192.168.x.x:8000` on a phone from your PC's development server generally will **not** grant camera access.

### Test on a desktop

From this folder, run `python -m http.server 8000`; open <http://localhost:8000> on the same PC. You need a webcam and to show the Hiro marker to it. `localhost` is treated as a potentially trustworthy origin for camera permissions. Double-clicking the HTML file using `file://` is unreliable.

## Replace with your own precomputed field lines

1. Export **visible tube/ribbon meshes** as GLB 2.0. `LineRenderer` objects or procedural lines in Unity often need to be converted/baked to meshes first; an exporter may otherwise omit them. Keep the field geometry compact (start with 10–30 selected lines).
2. Put the file in `models/` and change in `index.html`:
   `gltf-model="url(models/demo_fieldlines.glb)"`
   to your new path.
3. Ensure **model dimensions and marker coordinates match**. In this demo, marker plane is local XZ, Y is normal to the marker; the marker is about one *AR scene unit* wide by convention, and the sample coil has radius 0.55 marker widths and center height 0.60 marker widths.
4. Set the `position`, `rotation`, and `scale` of `<a-entity id="field-anchor">` to align your model with the apparatus. You can bake the full marker-to-coil rigid transformation into your GLB instead.
5. Validate with a simple calibration overlay first: draw your reference loop in AR and check whether its center, diameter and tilt match the actual coil from several viewpoints.

**Example calibration:** if the black square measures 10 cm and the coil center is 8 cm to the right of the marker center in local +X, the approximate displacement is `position="0.8 0 0"`, *provided your GLB geometry is already in marker-width units*. If the model uses **meters**, normalize its geometry by the marker's physical width in meters, e.g. a 0.08 m offset corresponds to 0.8 normalized units for a 0.10 m marker. Confirm axes and scale on-screen with a ruler; the `size` and camera intrinsics can affect reported registration.

## Important pitfalls

- **The overlay is simulated, not sensed.** The phone does not measure the coil's magnetic field; the drawing is based on supplied model and current.
- **Marker should stay visible** for robust tracking; occlusion by a hand or the equipment may make the field disappear. Marker-tracking jitter is normal.
- **The marker must be mechanically fixed relative to the coil.** A loose marker on the table is insufficient for quantitative comparisons.
- **No actual occlusion of virtual lines by the real coil** is implemented. Lines can be drawn in front of the physical conductor even where they should be behind it.
- **GLB line primitives may render as thin 1-pixel lines.** Tube meshes provide more consistent line thickness on mobile.
- **Camera and model loading are browser-dependent.** Test on an Android Chrome phone and an iPhone Safari before classroom use.
- **Static hosting:** assets should share the same host to avoid cross-origin problems. The external JavaScript libraries and the marker image require Internet; bundle pinned library scripts locally for fully offline operation.

## Future enhancements

- Add switchable field-line / component-vector / point-probe representations.
- Link current `I` to line thickness or displayed numerical magnitude; the **field-line geometry** of a linear current-loop model does not depend on the value of `I` (when the current direction is unchanged).
- Provide a model-only 3D viewer for unsupported devices or students who do not wish to use the camera.
- Consider a **custom AR.js pattern marker** if you want to include apparatus identification and calibration references.

## Sources

- AR.js marker-based guide: <https://ar-js-org.github.io/AR.js-Docs/marker-based/>
- AR.js marker-based example: <https://github.com/AR-js-org/AR.js/blob/master/aframe/examples/marker-based/basic.html>
- A-Frame GLB model component: <https://aframe.io/docs/1.8.0/components/gltf-model.html>
- Camera HTTPS requirement: <https://developer.mozilla.org/en-US/docs/Web/API/MediaDevices/getUserMedia>
- GitHub Pages setup: <https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site>
"# field-ar" 
