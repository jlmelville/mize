---
title: "Quasi-Hyperbolic Momentum"
date: "April 30, 2020"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
editor_options: 
  markdown: 
    wrap: 72
---

*December 12 2021*: reworded all uses of "generalized" momentum (a term of my
own bad invention) with "unified momentum" as used by [Zou and
co-workers](https://arxiv.org/abs/1808.03408).

When writing about [Nesterov Accelerated Gradient and
Momentum](https://jlmelville.github.io/mize/nesterov.html) I wrote out a
formula that made it easy to interpolate between Nesterov accelerated
gradient (NAG) and classical momentum (CM), and then discovered that
[Quasi-Hyperbolic Momentum](https://arxiv.org/abs/1810.06801) seemed to
do the same thing. It is not immediately clear from the presentation in
the QHM paper if the expressions are exactly the same, so here are some
notes on my efforts to rewrite QHM into a form closer to the NAG page.

A brief definition of the symbols I use:

-   $\theta$: the parameters being optimized.
-   $\varepsilon$: the learning rate or step size.
-   $f$: the cost or loss function to be minimized.
-   $\nabla f\left(\theta \right)$: the gradient of the cost function
    with respect to the parameters.
-   $\mu$: the momentum coefficient.
-   $v$: the velocity or update vector.
-   $s$: the steepest descent step.

Any of these can vary with iteration which is marked by the subscript
$t$. Most importantly for this discussion, the update vector is always:

$$
v_t \overset{\mathrm{def}}{=} \theta_t - \theta_{t-1}
$$

it can be expressed in other ways depending on the type of momentum
being used, but that definition always holds. Also, the steepest descent
step is:

$$
s_t \overset{\mathrm{def}}{=} -\varepsilon_t \nabla f\left(\theta_t \right)
$$

Some others will be introduced as I need them. I used equivalence symbol
($\equiv$) for definitions on the NAG page, but I thought I would change
it up and use $\overset{\mathrm{def}}{=}$ for this page, just to see
which one looks nicer on the screen. I'm sure I am using them both
incorrectly, but hopefully the intention comes across.

## Standard Definitions

### Classical Momentum

$$
\theta_{t+1} = \theta_t + s_t + \mu v_t
$$

### Nesterov Accelerated Gradient

$$
\theta_{t+1} = \theta_t + \left(1 + \mu \right) s_t + \mu v_t - \mu s_{t-1}
$$

### Unified Momentum

From the above, it's pretty easy to come up with a scheme to interpolate
between CM and NAG:

$$
\theta_{t+1} = \theta_t + \left(1 +  \beta \mu \right) s_t + \mu v_t - \beta \mu s_{t-1}
$$

where $\beta = 0$ gives CM and $\beta = 1$ gives NAG. The earliest publication
I am aware of that gives this expression is [Zou and
co-workers](https://arxiv.org/abs/1808.03408) from 2018. I'd love to know of any
earlier uses.

## Damped Momentum

In the [QHM](https://arxiv.org/abs/1810.06801) paper, a "damped" or
"normalized" expression for momentum is given. This reweights the
expression so that the weighting coefficients of the gradient descent
and the momentum term sums to one. The idea here is to stop the length
of the update, $v_t$ from being dependent on the momentum, i.e. in the
above expressions, the more momentum you have, the longer the step sizes
will be. This confounds any analysis of how momentum helps: better
direction or just a longer step size?

Damping CM just means weighting the gradient descent by $(1 - \mu)$. For
NAG, this requires weighing $s_t$ (and $s_{t-1}$) by $(1 - \mu)$ in all
parts of the expression to take into account that it is making a smaller
contribution everywhere.

### Damped Classical Momentum

$$
\theta_{t+1} = \theta_t + \left(1 - \mu \right) s_t + \mu v_t
$$

### Damped Nesterov Accelerated Gradient

$$
\theta_{t+1} = \theta_t + \left(1 - \mu^2 \right) s_t 
+ \mu v_t 
- \left(1 - \mu \right) \mu s_{t-1}
$$

or:

$$
\theta_{t+1} = \theta_t + \left(1 + \mu \right) \left(1 - \mu \right) s_t 
+ \mu v_t 
- \left(1 - \mu \right) \mu s_{t-1}
$$

### Damped Unified Momentum

$$
\theta_{t+1} = \theta_t + 
\left( 1 + \beta \mu \right)\left(1 - \mu \right) s_t 
+ \mu v_t 
- \beta \left(1 - \mu \right) \mu s_{t-1}
$$

The advantage of the damped expressions is that if you "unroll" an
expression (with gradient descent on the first step) the update vector
is just a weighted average of previous updates, i.e. the weighting
factors sum to one without any normalization needed.

### Quasi-Hyperbolic Momentum

$$
\theta_{t+1} = \theta_t
+ \left( 1 - \nu \mu \right ) s_t
+ \mu v_t
- \left( 1 - \nu \right ) \mu s_{t-1}
$$

$\nu = 0$ gives gradient descent, $\nu = 1$ gives damped classical
momentum, and $\nu = \mu$ gives damped NAG.

### Deriving the QHM update

The QHM expression above seems to differ from the one given in the [QHM
paper](https://arxiv.org/abs/1810.06801). Also, it's not immediately
apparent that $\nu = 0$ gives gradient descent, despite my assertion.

Here we'll start with QHM update in a way that most resembles how the
paper expresses it and then rewrite it into the form above. The QHM
update is in a form that is both "damped" *and* assumes a constant
learning rate, that latter allowing the extraction of the $\varepsilon$
out of the various update rules. Don't worry about why the QHM authors
want to do that for now.

Also, I need to introduce a couple of new symbols.

First, we won't be evaluating the gradient of the objective function $f$
anywhere except at $\theta_t$, so let's define a short-hand to represent
the gradient at iteration $t$. It makes perfect sense to me to represent
the gradient using the symbol $g$, but the QHM authors had other plans,
and use that to represent the "momentum buffer" (see below). To avoid
confusion with the notation in the paper I will reluctantly choose a
different symbol. The QHM authors use
$\nabla \hat{L} \left( \theta_t \right)$ for the gradient of the
objective function, so let's define:

$$
l_t \overset{\mathrm{def}}{=} \nabla f\left(\theta_t \right)
$$

Second, QHM uses the expression close to the (damped) classical momentum
update, which is calls the "momentum buffer". But note that it is *not*
the update vector used in QHM, so we can't call it $v_t$. The QHM
authors refer to this as $g_t$, but I refuse on principle to use this
notation as I am incapable of not mistaking that for the gradient
whenever I read it. Instead we'll call it $m_t$ (as in m for momentum):

$$
m_{t+1} = (1 - \mu) l_t + \mu m_t
$$

We are now ready to see the QHM update:

$$
\theta_{t+1} = \theta_t 
- \varepsilon 
\left[
\left( 1 - \nu \right ) l_t 
+ \nu m_{t+1}
\right]
$$

When $\nu = 0$, we get back plain old gradient descent. That is quite
easy to see. What's less easy to see is that $\nu = \mu$ gets you damped
NAG and $\nu = 1$ gets you damped classical momentum. Now is the time to
do some jiggery-pokery to the QHM equation and see what comes out.

#### Rewriting QHM

The update vector for QHM is:

$$
v_{t+1} =
- \varepsilon 
\left[
\left( 1 - \nu \right ) l_t + \nu m_{t+1}
\right]
$$

Let's expand:

$$
v_{t+1} =
- \varepsilon \left( 1 - \nu \right ) l_t
- \varepsilon \nu m_{t+1}
$$

and rearrange for $m_{t+1}$:

$$
m_{t+1} =
- \frac{\left( 1 - \nu \right ) l_t}{\nu}
- \frac{v_{t+1}}{\varepsilon \nu}
$$

and then subtract one from each subscript to get an expression for
$m_t$:

$$
m_t =
- \frac{\left( 1 - \nu \right ) l_{t-1}}{\nu}
- \frac{v_t}{\varepsilon \nu}
$$

Next, let's substitute the $m_{t+1}$ in the QHM update with the momentum
buffer update expression:

$$
\theta_{t+1} = \theta_t 
- \varepsilon 
\left[
\left( 1 - \nu \right ) l_t 
+ \nu (1 - \mu) l_t 
+ \nu \mu m_t
\right]
$$

and group the $l_t$ terms together to get:

$$
\theta_{t+1} = \theta_t 
- \varepsilon 
\left[
\left( 1 - \nu \mu \right ) l_t 
+ \nu \mu m_t
\right]
$$

Now we can substitute in the expression for $m_t$ we just came up with:

$$
\theta_{t+1} = \theta_t 
- \varepsilon 
\bigg\{
\left( 1 - \nu \mu \right ) l_t 
+ \nu \mu 
\left[
- \frac{\left( 1 - \nu \right ) l_{t-1}}{\nu}
- \frac{v_t}{\varepsilon \nu}
\right]
\bigg\}
$$

to give:

$$
\theta_{t+1} = \theta_t 
- \varepsilon 
\left[
\left( 1 - \nu \mu \right ) l_t 
- \mu \left( 1 - \nu \right ) l_{t-1}
- \mu \frac{v_t}{\varepsilon}
\right]
$$

Finally, let's expand the bracket, using $s_t = -\varepsilon l_t$ to
give this three-term QHM expression:

$$
\theta_{t+1} = \theta_t 
+ \left( 1 - \nu \mu \right ) s_t
- \mu \left( 1 - \nu \right ) s_{t-1}
+ \mu v_t
$$

Hopefully it is now easier to see that if you subsitute $\nu = \mu$, you
get damped NAG, and $\nu = 1$ gets you damped classical momentum. With
$\nu = 0$, when QHM should give us gradient descent, the expression
looks a little odd as there are still terms that contain $\mu$:

$$
\theta_{t+1} = \theta_t 
+ s_t
- \mu s_{t-1}
+ \mu v_t
$$

If you wade through the expressions above, you'll see that $v_t$ can be
written as:

$$
v_t =
\left( 1 - \nu \right ) s_{t-1}
- \varepsilon \nu m_t
$$

And when $\nu = 0$ the second term disappears so we get:

$$
\mu v_t = \mu s_{t-1}
$$

Therfore, for $\nu = 0$, the second and third terms in the three-term
QHM equation cancel, leaving just gradient descent as expected.

#### Undamped QHM?

I have been unsuccessful in coming up with an undamped version of QHM
which generalizes to NAG, CM and steepest descent. Which is not to say
it can't be done.

### Constant Learning Rate

QHM goes one step further to also assume a fixed learning rate
$\varepsilon$ at each step. The advantage of *that* can be seen by
unrolling the classical momentum (undamped to be as clear as possible)
for a few steps:

$$
\theta_2 = \theta_1 + s_1 \\
\theta_2 = \theta_1 - \varepsilon l_1 \\
\implies v_2 = - \varepsilon l_1
$$

$$
\theta_3 = \theta_2 + s_2 + \mu v_2 \\
\theta_3 = 
\theta_2 - \varepsilon l_2
- \mu \varepsilon l_1 \\
\theta_3 =
\theta_2 - \varepsilon \left( l_2 + \mu l_1 \right ) \\
\implies v_3 = - \varepsilon \left( l_2 + \mu l_1 \right)
$$

$$
\theta_4 = \theta_3 + s_3 + \mu v_3
\\
\theta_4 = \theta_3 - \varepsilon l_3
- \mu \left[ - \varepsilon \left( l_2 + \mu l_1 \right) \right]
\\
\implies v_4 = - \varepsilon \left[ l_3 + \mu \left( l_2 + \mu l_1 \right) \right]
$$

You can keep writing out the updates, but with a constant $\varepsilon$
you can always rearrange so that the update is in the form
$v_t = -\varepsilon d_t$ where $d_t$ is just the direction of the
update. This is also true for NAG and the damped versions. The
motivation to do this for QHM is to write the update as a weighted
average of gradients so the update can also be seen as a gradient
estimator in a stochastic gradient descent context.

The other advantage of this is that in terms of translating from the
un-damped expressions to the damped forms, we can always start with the
un-damped expression using a momentum coefficient $\mu^*$ and then
transform into the damped expression by rescaling the contributions from
the gradient descent and momentum, and then rescaling the learning rate,
to compensate, i.e.:

$$
\theta_{t+1} = \theta_t -\varepsilon \left(l_t + \mu^* v_t \right) = \\
\theta_t -\frac{\varepsilon}{1 + \mu^*} \left[ \left(1 - \mu \right) l_t + \mu v_t \right] = \\
\theta_t - \left(1 - \mu \right) \varepsilon \left[ \left(1 - \mu \right) l_t + \mu v_t \right]
$$

with

$$
\mu = \frac{\mu^*}{1 + \mu^*}
$$

As noted in the paper, this shrinks the update by a factor of $1 - \mu$
and this can be compensated by increasing $\varepsilon$ accordingly.

Note also that the recommended setting of the damped momentum
coefficient $\mu = 0.999$ in QHM would therefore translate to a very
large undamped momentum coefficient of $\mu^* = 999$.

### Unified Momentum Again

For completeness (and my own edification), we can take the same approach
with unified momentum and have the learning rate apply to the entire
direction of the update ($d_t$) rather than just the steepest descent
direction and therefore get an update that is a weighted average of
gradients. To make life easier notationally by removing lots of negative
signs, here's one more definition, the steepest descent direction:
$p_t \overset{\mathrm{def}}{=} -\nabla f\left(\theta_t \right)$.

The resulting unified momentum expression below also allows for
varying pretty much everything at each step: $\mu$, $\varepsilon$ and
$\beta$ although I don't have any suggestions for how you would go about
tuning these values:

$$
\theta_{t+1} = \theta_t + \varepsilon_t \left[
\left( 1 + \beta_t \mu_t \right)\left(1 - \mu_t \right) p_t 
+ \mu_t d_{t-1} 
- \beta_{t-1} \left(1 - \mu_{t-1} \right) \mu_t p_{t-1}
\right]
$$

## Sort of Related Ideas

Papers that (mainly) also focus on momentum.

-   [Understanding the Role of Momentum in Stochastic Gradient
    Methods](http://papers.nips.cc/paper/9158-understanding-the-role-of-momentum-in-stochastic-gradient-methods)

Lots more on the theoretical performance of QHM and a wealth of
supplementary material.

-   [Gradient descent with momentum --- to accelerate or to
    super-accelerate?](https://arxiv.org/abs/2001.06472)

Uses the Sutskever expression for NAG but ups the "look-ahead" concept
by taking *two* steps. The first step is the momentum step, then a
second step is taken using a large value of $\mu$ (the authors recommend
$\mu = 5$ at the start of the optimization, and then $\mu = 2$ later),
and the gradient is evaluated at that latter point. The gradient descent
step still takes place at the end point of the original moment step, but
uses the gradient evaluated at the second location.

-   [Decaying momentum helps neural network
    training](https://arxiv.org/abs/1910.04952)

This proposes a decaying momentum ("Demon") schedule of:

$$
\mu_t = \mu_0  \frac{1 - t / T}{\left( 1 - \mu_0 \right) + \mu_0 \left( 1 - t/T \right)}
$$

where $\mu_0$ is the initial momentum and $T$ is the total number of
iterations.

-   [Negative Momentum for Improved Game
    Dynamics](https://arxiv.org/abs/1807.04740)

This paper shows a trend in recent GAN papers of increasingly small
value for momentum and then proposes negative momentum. The analogy here
is the idea of friction in damping down oscillations. No connection with
anything discussed here, just an example of how weird momentum is.

-   [An updated overview of recent gradient descent
    algorithms](https://johnchenresearch.github.io/demon/)

Details and benchmarks a lot of the methods and papers listed here as
well as several other adaptive learning rate methods (e.g. Adam and its
many, many variants).

-   [YellowFin and the Art of Momentum
    Tuning](https://arxiv.org/abs/1706.03471)

Describes an adaptive momentum (and learning rate) system.

-   [Online Learning Rate Adaptation with Hypergradient
    Descent](https://arxiv.org/abs/1703.04782)

This has nothing to do with momentum but is another approach to
automatically choosing a learning rate. Perhaps its approach will be
applied to momentum tuning one day.
