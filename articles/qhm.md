# Quasi-Hyperbolic Momentum

*December 12 2021*: reworded all uses of “generalized” momentum (a term
of my own bad invention) with “unified momentum” as used by [Zou and
co-workers](https://arxiv.org/abs/1808.03408). Also, pointed to the
expression for Nesterov momentum in the [Nadam
paper](https://openreview.net/forum?id=OM0jvwB8jIp57ZJjtNEZ).

When writing about [Nesterov Accelerated Gradient and
Momentum](https://jlmelville.github.io/mize/articles/nesterov.html) I
wrote out a formula that made it easy to interpolate between Nesterov
accelerated gradient (NAG) and classical momentum (CM), and then
discovered that [Quasi-Hyperbolic
Momentum](https://arxiv.org/abs/1810.06801) seemed to do the same thing.
It is not immediately clear from the presentation in the QHM paper if
the expressions are exactly the same, so here are some notes on my
efforts to rewrite QHM into a form closer to the NAG page.

A brief definition of the symbols I use:

- \theta: the parameters being optimized.
- \varepsilon: the learning rate or step size.
- f: the cost or loss function to be minimized.
- \nabla f(\theta): the gradient of the cost function with respect to
  the parameters.
- \mu: the momentum coefficient.
- v: the velocity or update vector.
- s: the steepest descent step.

Any of these can vary with iteration which is marked by the subscript t.
Most importantly for this discussion, the update vector is always:

v_t \overset{\mathrm{def}}{=} \theta_t - \theta\_{t-1}

it can be expressed in other ways depending on the type of momentum
being used, but that definition always holds. Also, the steepest descent
step is:

s_t \overset{\mathrm{def}}{=} -\varepsilon_t \nabla f(\theta_t)

## Standard Definitions

### Classical Momentum

\theta\_{t+1} = \theta_t + s_t + \mu v_t

### Nesterov Accelerated Gradient

\theta\_{t+1} = \theta_t + (1 + \mu) s_t + \mu v_t - \mu s\_{t-1}

### Unified Momentum

From the above, it’s pretty easy to come up with a scheme to interpolate
between CM and NAG:

\theta\_{t+1} = \theta_t + (1 + \beta \mu) s_t + \mu v_t - \beta \mu
s\_{t-1}

where \beta = 0 gives CM and \beta = 1 gives NAG. The earliest
publication I am aware of that gives this expression is [Zou and
co-workers](https://arxiv.org/abs/1808.03408) from 2018. I’d love to
know of any earlier uses.

## Damped Momentum

In the [QHM](https://arxiv.org/abs/1810.06801) paper, a “damped” or
“normalized” expression for momentum is given. This reweights the
expression so that the weighting coefficients of the gradient descent
and the momentum term sum to one. The idea here is to stop the length of
the update, v_t from being dependent on the momentum, i.e. in the above
expressions, the more momentum you have, the longer the step sizes will
be. This confounds any analysis of how momentum helps: better direction
or just a longer step size?

Damping CM just means weighting the gradient descent by (1 - \mu). For
NAG, this requires weighting s_t (and s\_{t-1}) by (1 - \mu) in all
parts of the expression to take into account that it is making a smaller
contribution everywhere.

### Damped Classical Momentum

\theta\_{t+1} = \theta_t + (1 - \mu) s_t + \mu v_t

### Damped Nesterov Accelerated Gradient

\theta\_{t+1} = \theta_t + (1 - \mu^2) s_t + \mu v_t - (1 - \mu) \mu
s\_{t-1}

or:

\theta\_{t+1} = \theta_t + (1 + \mu) (1 - \mu) s_t + \mu v_t - (1 - \mu)
\mu s\_{t-1}

### Damped Unified Momentum

\theta\_{t+1} = \theta_t + (1 + \beta \mu)(1 - \mu) s_t + \mu v_t -
\beta (1 - \mu) \mu s\_{t-1}

The advantage of the damped expressions is that if you “unroll” an
expression (with gradient descent on the first step) the update vector
is just a weighted average of previous updates, i.e. the weighting
factors sum to one without any normalization needed.

### Quasi-Hyperbolic Momentum

\theta\_{t+1} = \theta_t + (1 - \nu \mu) s_t + \mu v_t - (1 - \nu) \mu
s\_{t-1}

\nu = 0 gives gradient descent, \nu = 1 gives damped classical momentum,
and \nu = \mu gives damped NAG.

### Deriving the QHM update

The QHM expression above seems to differ from the one given in the [QHM
paper](https://arxiv.org/abs/1810.06801). Also, it’s not immediately
apparent that \nu = 0 gives gradient descent, despite my assertion.

Here we’ll start with QHM update in a way that most resembles how the
paper expresses it and then rewrite it into the form above. The QHM
update is in a form that is both “damped” *and* assumes a constant
learning rate, that latter allowing the extraction of the \varepsilon
out of the various update rules. Don’t worry about why the QHM authors
want to do that for now.

I’ll keep \mu and \nu fixed too, and take \varepsilon \> 0 for this
derivation.

Also, I need to introduce a couple of new symbols.

First, we won’t be evaluating the gradient of the objective function f
anywhere except at \theta_t, so let’s define a short-hand to represent
the gradient at iteration t. It makes perfect sense to me to represent
the gradient using the symbol g, but the QHM authors had other plans,
and use that to represent the “momentum buffer” (see below). To avoid
confusion with the notation in the paper I will reluctantly choose a
different symbol. The QHM authors use \nabla \hat{L} (\theta_t) for the
gradient of the objective function, so let’s define:

\ell_t \overset{\mathrm{def}}{=} \nabla f(\theta_t)

Second, QHM uses the expression close to the (damped) classical momentum
update, which it calls the “momentum buffer”. But note that it is *not*
the update vector used in QHM, so we can’t call it v_t. The QHM authors
refer to this as g_t, but I refuse on principle to use this notation as
I am incapable of not mistaking that for the gradient whenever I read
it. Instead we’ll call it m_t (as in m for momentum):

m\_{t+1} = (1 - \mu) \ell_t + \mu m_t

We are now ready to see the QHM update:

\theta\_{t+1} = \theta_t - \varepsilon \left\[ (1 - \nu) \ell_t + \nu
m\_{t+1} \right\]

When \nu = 0, we get back plain old gradient descent. That is quite easy
to see.

*December 12 2021*: Before contorting the QHM equation further, the rest
of this section is new, mainly to draw attention to how QHM also
generates Dozat’s version of NAG in the Nadam paper.

With \nu = 1, the only thing left in the square brackets is the momentum
buffer m\_{t+1} which combines with \varepsilon to give the damped
classical momentum update. That’s also quite easy to see.

What’s less easy to see is that \nu = \mu gets you damped NAG. I’d like
to draw your attention to the immediate result:

\theta\_{t+1} = \theta_t - \varepsilon \left\[ (1 - \mu) \ell_t + \mu
m\_{t+1} \right\]

The bit in the square brackets has an identical form to the momentum
buffer except that m\_{t} has been replaced by m\_{t+1}. For some reason
I find this delightful:

- You want damped classical momentum? Form a weighted average of the
  current gradient with the previous iteration’s damped classical
  momentum direction.
- You wanted damped NAG? Just do the above, then do it *again*: form a
  weighted average of the current gradient with the new damped classical
  momentum direction from the step above.

As I have apparently made a hobby of collecting different researchers
giving different expressions for Nesterov momentum, here I must point
out that this is the form of the Nesterov momentum expression given by
[Dozat in 2016](https://openreview.net/forum?id=OM0jvwB8jIp57ZJjtNEZ):
see equation 11.

Anyway, substituting in m\_{t+1}, we get to:

\theta\_{t+1} = \theta_t - \varepsilon \left\[ (1 + \mu) (1 - \mu)
\ell_t + \mu^2 m_t \right\]

That’s closer to the damped NAG result, but there is a conspicuous
absence of contributions from the gradient step at the previous
iteration, and it’s still written in terms of the momentum buffer and
not the NAG update vector (although that’s not necessarily a bad thing).
Now is the time to do some jiggery-pokery to the QHM equation and see
what comes out.

#### Rewriting QHM

The update vector for QHM is:

v\_{t+1} = - \varepsilon \left\[ (1 - \nu) \ell_t + \nu m\_{t+1}
\right\]

Let’s expand:

v\_{t+1} = - \varepsilon (1 - \nu) \ell_t - \varepsilon \nu m\_{t+1}

and, for \nu \ne 0, rearrange for m\_{t+1}:

m\_{t+1} = - \frac{(1 - \nu) \ell_t}{\nu} - \frac{v\_{t+1}}{\varepsilon
\nu}

and then subtract one from each subscript to get an expression for m_t:

m_t = - \frac{(1 - \nu) \ell\_{t-1}}{\nu} - \frac{v_t}{\varepsilon \nu}

Next, let’s substitute the m\_{t+1} in the QHM update with the momentum
buffer update expression:

\theta\_{t+1} = \theta_t - \varepsilon \left\[ (1 - \nu) \ell_t + \nu
(1 - \mu) \ell_t + \nu \mu m_t \right\]

and group the \ell_t terms together to get:

\theta\_{t+1} = \theta_t - \varepsilon \left\[ (1 - \nu \mu) \ell_t +
\nu \mu m_t \right\]

Now we can substitute in the expression for m_t we just came up with:

\theta\_{t+1} = \theta_t - \varepsilon \bigg\\ (1 - \nu \mu) \ell_t +
\nu \mu \left\[ - \frac{(1 - \nu) \ell\_{t-1}}{\nu} -
\frac{v_t}{\varepsilon \nu} \right\] \bigg\\

to give:

\theta\_{t+1} = \theta_t - \varepsilon \left\[ (1 - \nu \mu) \ell_t -
\mu (1 - \nu) \ell\_{t-1} - \mu \frac{v_t}{\varepsilon} \right\]

Finally, let’s expand the bracket, using s_t = -\varepsilon \ell_t to
give this three-term QHM expression:

\theta\_{t+1} = \theta_t + (1 - \nu \mu) s_t - \mu (1 - \nu) s\_{t-1} +
\mu v_t

Hopefully it is now easier to see that if you substitute \nu = \mu, you
get damped NAG, and \nu = 1 gets you damped classical momentum. We
divided by \nu on the way here, so let’s check \nu = 0 directly.
Substituting it into the three-term expression leaves something that
looks a little odd, as there are still terms that contain \mu:

\theta\_{t+1} = \theta_t + s_t - \mu s\_{t-1} + \mu v_t

For comparison, the original QHM update at the previous iteration gives:

v_t = (1 - \nu) s\_{t-1} - \varepsilon \nu m_t

And when \nu = 0 the second term disappears so we get:

\mu v_t = \mu s\_{t-1}

Therefore, for \nu = 0, the second and third terms in the three-term QHM
equation cancel, leaving just gradient descent as expected.

#### Undamped QHM?

I have been unsuccessful in coming up with an undamped version of QHM
which generalizes to NAG, CM and steepest descent. Which is not to say
it can’t be done.

### Constant Learning Rate

QHM goes one step further to also assume a fixed learning rate
\varepsilon at each step. The advantage of *that* can be seen by
unrolling the classical momentum (undamped to be as clear as possible)
for a few steps:

\begin{aligned} \theta_2 &= \theta_1 + s_1 \\ &= \theta_1 - \varepsilon
\ell_1 \\ v_2 &= -\varepsilon \ell_1 \end{aligned}

\begin{aligned} \theta_3 &= \theta_2 + s_2 + \mu v_2 \\ &= \theta_2 -
\varepsilon \ell_2 - \mu \varepsilon \ell_1 \\ &= \theta_2 - \varepsilon
(\ell_2 + \mu \ell_1) \\ v_3 &= -\varepsilon (\ell_2 + \mu \ell_1)
\end{aligned}

\begin{aligned} \theta_4 &= \theta_3 + s_3 + \mu v_3 \\ &= \theta_3 -
\varepsilon \ell_3 + \mu \left\[-\varepsilon (\ell_2 + \mu
\ell_1)\right\] \\ v_4 &= -\varepsilon \left\[\ell_3 + \mu (\ell_2 + \mu
\ell_1)\right\] \end{aligned}

You can keep writing out the updates, but with a constant \varepsilon
you can always rearrange so that the update is in the form v\_{t+1} =
\varepsilon d_t where d_t is just the direction of the update. This is
also true for NAG and the damped versions. The motivation to do this for
QHM is to express its square-bracketed term as a weighted average of
gradients, which can also be seen as a gradient estimator in a
stochastic gradient descent context.

We can also make the connection between damped and undamped classical
momentum explicit. Keep 0 \leq \mu \< 1 fixed, and use u_t for an
undamped buffer of gradients:

u\_{t+1} = \ell_t + \mu u_t

The damped buffer m_t introduced above uses the same momentum
coefficient:

m\_{t+1} = (1 - \mu) \ell_t + \mu m_t

Starting both buffers at zero gives m_t = (1 - \mu) u_t. So we can write
the classical momentum update in either form:

\begin{aligned} \theta\_{t+1} &= \theta_t -
\varepsilon\_{\mathrm{undamped}} u\_{t+1} \\ &= \theta_t -
\frac{\varepsilon\_{\mathrm{undamped}}}{1 - \mu} m\_{t+1} \end{aligned}

As noted in the paper, damping shrinks the update by a factor of 1 - \mu
if we keep the learning rate the same. To recover the undamped updates,
choose:

\varepsilon\_{\mathrm{damped}} =
\frac{\varepsilon\_{\mathrm{undamped}}}{1 - \mu}

The momentum coefficient stays at \mu. For example, with \mu = 0.999,
the damped learning rate needs to be 1000 times the undamped learning
rate to give the same classical momentum updates.

### Unified Momentum Again

For completeness (and my own edification), we can take the same approach
with unified momentum and have the learning rate apply to the entire
direction of the update (d_t) rather than just the steepest descent
direction and therefore get an update that is a weighted average of
gradients. To make life easier notationally by removing lots of negative
signs, here’s one more definition, the steepest descent direction: p_t
\overset{\mathrm{def}}{=} -\nabla f(\theta_t).

The resulting unified momentum expression below also allows for varying
pretty much everything at each step: \mu, \varepsilon and \beta although
I don’t have any suggestions for how you would go about tuning these
values:

\theta\_{t+1} = \theta_t + \varepsilon_t \left\[ (1 + \beta_t \mu_t)(1 -
\mu_t) p_t + \mu_t d\_{t-1} - \beta\_{t-1} (1 - \mu\_{t-1}) \mu_t
p\_{t-1} \right\]

## Sort of Related Ideas

Papers that (mainly) also focus on momentum.

- [Understanding the Role of Momentum in Stochastic Gradient
  Methods](https://papers.nips.cc/paper/9158-understanding-the-role-of-momentum-in-stochastic-gradient-methods)

Lots more on the theoretical performance of QHM and a wealth of
supplementary material.

- [Gradient descent with momentum — to accelerate or to
  super-accelerate?](https://arxiv.org/abs/2001.06472)

Uses the Sutskever expression for NAG but ups the “look-ahead” concept
by taking *two* steps. The first step is the momentum step, then a
second step is taken using a large value of \mu (the authors recommend
\mu = 5 at the start of the optimization, and then \mu = 2 later), and
the gradient is evaluated at that latter point. The gradient descent
step still takes place at the end point of the original moment step, but
uses the gradient evaluated at the second location.

- [Decaying momentum helps neural network
  training](https://arxiv.org/abs/1910.04952)

This proposes a decaying momentum (“Demon”) schedule of:

\mu_t = \mu_0 \frac{1 - t/T}{(1 - \mu_0) + \mu_0(1 - t/T)}

where \mu_0 is the initial momentum and T is the total number of
iterations.

- [Negative Momentum for Improved Game
  Dynamics](https://arxiv.org/abs/1807.04740)

This paper shows a trend in recent GAN papers of increasingly small
value for momentum and then proposes negative momentum. The analogy here
is the idea of friction in damping down oscillations. No connection with
anything discussed here, just an example of how weird momentum is.

- [An updated overview of recent gradient descent
  algorithms](https://johnchenresearch.github.io/demon/)

Details and benchmarks a lot of the methods and papers listed here as
well as several other adaptive learning rate methods (e.g. Adam and its
many, many variants).

- [YellowFin and the Art of Momentum
  Tuning](https://arxiv.org/abs/1706.03471)

Describes an adaptive momentum (and learning rate) system.

- [Online Learning Rate Adaptation with Hypergradient
  Descent](https://arxiv.org/abs/1703.04782)

This has nothing to do with momentum but is another approach to
automatically choosing a learning rate. Perhaps its approach will be
applied to momentum tuning one day.
