# Visualizing networks with cograph

## Drawing a fitted network

A fitted `psychnets` network is a weighted graph of nodes and signed
edges. The picture of that graph is where its structure becomes
readable: which nodes connect, how strongly, and with what sign. Every
estimator in `psychnets` returns an object of class
`c("psychnet", "cograph_network")`, so a fitted network moves into the
`cograph` drawing engine with no conversion. The same model is then
rendered with different layouts, node aesthetics, edge styling,
thresholds, themes, and export settings, all from one object.

This vignette shows the plotting workflow from the default graph to a
figure ready for a paper.
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
is the drawing verb throughout. The chunks that call it are guarded, so
they run only when `cograph` is installed.

## The data and the fit

The bundled `SRL_GPT` data hold 300 learners scored on five
self-regulated -learning constructs (CSU, IV, SE, SR, TA).
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
estimates a Gaussian graphical model whose edges are partial
correlations selected by the extended Bayesian information criterion.

``` r

fit <- ebic_glasso(SRL_GPT)
fit
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 2.21e-10
```

## The default graph

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a `psychnet`
object delegates to
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
so a bare `plot(fit)` draws the network. Each node is a construct and
each edge is a partial correlation between two constructs given the
other three. The colour of an edge encodes the sign of that partial
correlation and the width encodes its magnitude, so a wide edge of one
colour is a strong positive association and a wide edge of the other
colour is a strong negative one.

``` r

plot(fit)
```

![](visualizing-with-cograph_files/figure-html/default-1.png)

Calling
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
directly draws the same graph and opens its full argument surface for
the customization below.

``` r

cograph::splot(fit)
```

![](visualizing-with-cograph_files/figure-html/splot-default-1.png)

## Layouts

The `layout` argument sets the algorithm that places the nodes. A
circular layout fixes the nodes on a ring, which keeps positions
comparable between figures; a spring layout places connected nodes near
one another, which makes clusters visible. The `seed` argument fixes the
random start of a stochastic layout so the figure is reproducible.

``` r

op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
cograph::splot(fit, layout = "circle", title = "circle")
cograph::splot(fit, layout = "spring", seed = 11, title = "spring")
```

![](visualizing-with-cograph_files/figure-html/layouts-1.png)

``` r

par(op)
```

## Nodes and labels

The node arguments set fill, border, label size, and shape. The
`scale_nodes_by` argument sizes each node by a centrality measure and
`node_size_range` sets the smallest and largest radius, so node area
reads as structural importance. Sizing the nodes by strength makes the
most connected construct the largest on the page.

``` r

cograph::splot(
  fit,
  layout = "circle",
  scale_nodes_by = "strength",
  node_size_range = c(4, 12),
  node_fill = c("#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"),
  node_border_color = "white",
  node_border_width = 2,
  label_size = 0.9,
  title = "Node size scaled by strength"
)
```

![](visualizing-with-cograph_files/figure-html/nodes-1.png)

## Edges

For a signed psychometric network, the edge encoding carries the
substantive result. The `edge_positive_color` and `edge_negative_color`
arguments set the two sign colours and `edge_width_range` sets the
mapping from magnitude to width. The `threshold` argument hides edges
below an absolute weight, which thins a dense graph down to its
strongest associations, and `edge_labels` prints the weight on each
retained edge.

``` r

cograph::splot(
  fit,
  layout = "circle",
  threshold = 0.05,
  edge_width_range = c(0.5, 5),
  edge_positive_color = "#2A9D8F",
  edge_negative_color = "#E76F51",
  edge_labels = TRUE,
  edge_label_size = 0.7,
  edge_label_bg = "white",
  title = "Thresholded edge weights"
)
```

![](visualizing-with-cograph_files/figure-html/edges-1.png)

[`cograph::plot_edge_weights()`](https://sonsoles.me/cograph/reference/plot_edge_weights.html)
draws the distribution of the edge weights on its own, which reads the
spread of associations and the location of the strong edges without the
graph layout.

``` r

cograph::plot_edge_weights(fit)
```

![](visualizing-with-cograph_files/figure-html/edge-weights-1.png)

## Bootstrap diagnostics

The bootstrap diagnostics are computed in `psychnets` and drawn by their
own base-R plot methods, covered fully in the companion vignette.
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
resamples the data and refits the network on each resample, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on its result
defaults to the edge-weight confidence intervals.

``` r

set.seed(1)
bs <- net_boot(SRL_GPT, method = "glasso", n_boot = 250, cores = 1)
plot(bs)
```

![](visualizing-with-cograph_files/figure-html/bootstrap-1.png)

The same bootstrap object holds the retained draws for the pairwise
edge-difference test. Under `type = "edge_diff"`,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
significance-box matrix of which edges differ from one another.

``` r

plot(bs, type = "edge_diff")
```

![](visualizing-with-cograph_files/figure-html/edge-diff-box-1.png)

[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
returns the pairwise differences as a tidy table, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) under
`style = "forest"` draws each difference as a point with its confidence
interval when the effect sizes matter more than the box display.

``` r

plot(difference_test(bs, type = "edge"), style = "forest")
```

![](visualizing-with-cograph_files/figure-html/edge-diff-forest-1.png)

## Themes

The `theme` argument applies a coordinated set of colour and styling
defaults without changing the fitted model, so the same network is
redrawn in a house style or a colourblind-safe palette from one call.

``` r

op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
cograph::splot(fit, theme = "minimal", title = "minimal")
cograph::splot(fit, theme = "colorblind", title = "colorblind")
```

![](visualizing-with-cograph_files/figure-html/themes-1.png)

``` r

par(op)
```

## Group networks

Passing `group =` to
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
fits one network per level of a grouping column and returns the
collection.
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the collection as a grid that shares one layout, so a node sits in
the same position in every panel and the panels are read side by side.

``` r

group_fit <- psychnet(grouped_srl, group = "source", method = "glasso")
cograph::splot(group_fit, layout = "circle", psych_styling = TRUE)
```

![](visualizing-with-cograph_files/figure-html/groups-1.png)

## Export

The `filename`, `width`, `height`, and `res` arguments write the figure
to a file at a chosen size and resolution when it is ready to save.

``` r

cograph::splot(fit, layout = "circle", filename = "network.png",
               width = 8, height = 8, res = 300)
```

## The splot arguments at a glance

| Task | Arguments |
|----|----|
| Layout | `layout`, `seed`, `layout_scale`, `layout_margin` |
| Nodes | `node_size`, `scale_nodes_by`, `node_size_range`, `node_fill`, `node_shape` |
| Labels | `labels`, `label_size`, `label_color`, `label_position` |
| Edges | `edge_width_range`, `edge_positive_color`, `edge_negative_color`, `edge_alpha` |
| Filtering | `threshold`, `minimum`, `maximum` |
| Edge labels | `edge_labels`, `edge_label_size`, `edge_label_bg`, `edge_label_style` |
| Themes | `theme`, `background`, `title` |
| Export | `filename`, `width`, `height`, `res` |
