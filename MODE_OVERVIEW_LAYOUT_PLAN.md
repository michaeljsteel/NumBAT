# Mode Overview Layout Plan

## Purpose

Improve the space efficiency of the combined mode overview while preserving
equal physical scaling in all field panels. The layout must work for both
near-square and highly elongated waveguides without changing requested field
content or individual-component plots.

## Current Layouts

The combined overview currently uses:

- a tall six-panel layout: two rows by three columns;
- a wide six-panel layout: three rows by two columns;
- a nine-panel layout when secondary components are displayed.

The two six-panel layouts have substantially different usable field-panel
aspect ratios. The existing ``wide_mode`` Boolean selects one manually, but
the default does not account for the visible waveguide aspect ratio.

## Layout Policy

Expose:

```python
overview_layout="auto"  # default
```

with explicit alternatives:

```python
overview_layout="tall"
overview_layout="wide"
```

``wide_mode`` remains a backwards-compatible override:

- ``wide_mode=True`` selects ``"wide"``;
- ``wide_mode=False`` selects ``"tall"``;
- ``wide_mode=None`` permits ``overview_layout`` or automatic selection.

When both options are supplied, explicit ``wide_mode`` takes precedence for
compatibility.

## Phase 1: Adaptive Six-Panel Selection

**Status: Complete**

For overviews that suppress secondary components:

1. Calculate the visible waveguide aspect ratio from the displayed x and y
   spans, including configured plot trimming.
2. Estimate each current template's field-cell aspect ratio from the figure
   dimensions and grid shape.
3. Select the template whose aspect is closest on a logarithmic ratio scale.
4. Keep the existing figure sizes and panel contents unchanged.

The nine-panel automatic layout remains the current tall arrangement. Explicit
``"wide"`` still selects the existing wide nine-panel arrangement.

## Phase 2: Bounded Adaptive Figure Dimensions

After observing Phase 1 outputs, optionally adjust the figure width/height
ratio while preserving a bounded total figure area. This should improve
extreme waveguide aspect ratios without generating unreasonably large files.

## Phase 3: Named Nine-Panel Mosaics

Add dedicated tall, square, and wide named mosaics for overviews with primary
and secondary Cartesian components. Use ``GridSpec`` or ``subplot_mosaic`` so
the mode-data summary, vector view, magnitude/intensity view, and related
Cartesian components remain semantically grouped.

## Compatibility and Validation

- Preserve the existing tall and wide templates as explicit choices.
- Preserve equal physical scaling in field panels.
- Preserve existing filenames, field selections, labels, and colour maps.
- Test automatic wide and tall choices, explicit overrides, and unchanged
  nine-panel automatic behavior.
