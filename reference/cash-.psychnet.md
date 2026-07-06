# Back-compatible field access for a psychnet object

The canonical (`str`-visible) fields are the lean netobject set; this
method adds virtual aliases so older/external accessors keep working
without storing redundant fields.

## Usage

``` r
# S3 method for class 'psychnet'
x$name
```

## Arguments

- x:

  A `psychnet` object.

- name:

  Field name. Canonical fields plus the legacy aliases `graph` (=
  `weights`), `labels` (= `nodes$label`), `n_nodes`, `n_edges`, and
  `n_obs` (= `n`).

## Value

The requested field, or `NULL` if neither a canonical field nor a known
alias.
