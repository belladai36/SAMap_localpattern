# Mathematical specification and guarantees

## 1. Input graph and aggregation

Let `W` be a nonnegative cross-species cell-cell edge-weight matrix, such as one
returned by SAMap. LocalSAMap conditions on `W`; it does not claim a new
guarantee for the procedure that learned the graph. For cell sets `A` and `B`,
define the mean edge score

\[
S(A,B)=\frac{1}{|A||B|}\sum_{i\in A}\sum_{j\in B}W_{ij}.
\]

This aggregation is descriptive. Inferential validity comes from the resampling
scheme and its assumptions, not from this formula alone.

## 2. Local module score

For standardized pseudobulk vectors `x` and `y` over a homologous module, and
nonnegative confidence weights `w`, define

\[
C_w(x,y)=\frac{\sum_g w_g x_g y_g}
 {\sqrt{\sum_g w_gx_g^2}\sqrt{\sum_g w_gy_g^2}}.
\]

**Proposition 1 (boundedness).** If both weighted norms are nonzero, then
`-1 <= C_w <= 1`.

**Proof.** Set `u_g=sqrt(w_g)x_g` and `v_g=sqrt(w_g)y_g`. The score is
`<u,v>/(||u|| ||v||)`. Cauchy-Schwarz gives
`|<u,v>| <= ||u|| ||v||`. Division by the positive denominator proves the
claim. QED.

The signed gene contribution is `c_g=w_g x_g y_g`; contributions sum to the
numerator, providing an exact additive explanation of that numerator.

## 3. Permutation inference

Let `T_obs` be a mapping statistic. Generate `R` label permutations that are
valid under the null, preferably within donor, experiment, and species strata.
With `K = sum_r 1(T_r >= T_obs)`, report

\[
\hat p=(K+1)/(R+1).
\]

**Theorem 1 (randomization-test validity).** If the observed labeling and the
`R` permuted labelings are exchangeable under the null and ties are handled by
the displayed weak inequality, then `P(p_hat <= alpha) <= alpha` on the
attainable grid (or after standard randomized tie breaking).

**Proof sketch.** Under exchangeability, the rank of `T_obs` among the `R+1`
statistics is uniform, up to conservative ties. The event `p_hat <= alpha`
contains at most `floor(alpha(R+1))` of the `R+1` ranks. Its probability is
therefore at most `alpha`. QED.

Exchangeability fails if labels are permuted across incompatible donors,
batches, developmental stages, or anatomical sampling designs. The software
cannot repair an invalid permutation design.

## 4. Monte Carlo precision

Conditional on the observed data, let `p_*` be the probability that a random
valid permutation is at least as extreme as observed, and let `K/R` estimate
it from independent permutations.

**Theorem 2 (Hoeffding bound).** For every `epsilon > 0`,

\[
P(|K/R-p_*|\ge\epsilon)\le 2\exp(-2R\epsilon^2).
\]

**Proof.** The exceedance indicators are independent Bernoulli variables in
`[0,1]`; apply Hoeffding's inequality to their mean. QED.

Thus, with probability at least `1-delta`, simulation error is no larger than
`sqrt(log(2/delta)/(2R))`. This concerns Monte Carlo error only.

## 5. Bootstrap recovery

Let `Z_b` indicate whether a prespecified mapping is recovered in bootstrap
replicate `b`, and let `pi_hat = B^{-1} sum_b Z_b`.

**Theorem 3 (conditional concentration).** If bootstrap replicates are
conditionally independent and identically generated, then

\[
P(|\hat\pi-E[Z_b\mid D]|\ge\epsilon\mid D)
\le 2\exp(-2B\epsilon^2).
\]

**Proof.** Conditional on data `D`, apply Hoeffding's inequality to the
Bernoulli recovery indicators. QED.

This estimates algorithmic stability under the chosen resampling design, not
the probability that a biological homology statement is true. Donor-level
bootstrap is preferred when donors are the independent experimental units.

## 6. Multiple testing

For `m` hypotheses with ordered p-values `p_(1) <= ... <= p_(m)`, the
Benjamini-Hochberg rule rejects through the largest `k` satisfying
`p_(k) <= k alpha/m`.

**Guarantee.** BH controls FDR at at most `alpha` when null p-values are valid
and independent, and under standard positive regression dependence (PRDS).
Arbitrary dependence requires a more conservative procedure such as
Benjamini-Yekutieli or an explicitly justified hierarchical design.

## 7. Specificity and optional coordinate agreement

For normalized nonnegative scores `q_b` over `K` candidate target types,

\[
H(q)=-\sum_b q_b\log(q_b)/\log K \in [0,1].
\]

Low entropy is specific; high entropy is ambiguous. The best-vs-second-best
margin is reported alongside it. If positions `p_a,p_b` are defensibly scaled
to `[0,1]`, optional coordinate agreement is

\[
A(a,b)=\exp(-(p_a-p_b)^2/(2\sigma^2)) \in (0,1].
\]

Coordinate agreement is reported separately from transcriptional similarity
because coordinates may shift across systems or may not be commensurable.

## 8. Claims this theory does not support

No theorem here proves that a high score is evolutionary homology. Results can
still be distorted by pseudoreplication, batch effects, stress/cell-cycle
programs, wrong annotations, incomplete proteomes, paralog substitution, or
incorrect coordinates. Conclusions require replication, marker inspection,
homology-graph sensitivity, and biological validation.
