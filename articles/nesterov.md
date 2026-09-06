# Nesterov Accelerated Gradient and Momentum

*December 16 2021*: I have rewritten a large chunk of this document to
make the derivations a little clearer and added a brief section on the
definition of Nesterov momentum used in the NAdam paper.

A way to express Nesterov Accelerated Gradient (NAG) in terms of a
regular momentum update was noted by [Sutskever and
co-workers](https://proceedings.mlr.press/v28/sutskever13.html), and
perhaps more importantly, when it came to training neural networks, it
seemed to work better than classical momentum schemes. This was further
confirmed by [Bengio and co-workers](https://arxiv.org/abs/1212.0901),
who provided an alternative formulation that might be easier to
integrate into existing software.

I implemented NAG as part of [mize](https://github.com/jlmelville/mize),
but found it was a bit more difficult than I’d anticipated to test that
I’d got the correct results.

It seems I am not alone in being confused by the exact details of
implementing NAG as a momentum scheme (e.g. comments and issues in the
deep learning projects
[keras](https://github.com/fchollet/keras/issues/966) and
[mocha.jl](https://github.com/pluskid/Mocha.jl/pull/47), which I found
from random googling), so mainly for the benefit of future me, I am
going to derive the Sutskever and Bengio formulations from NAG in
tedious detail. I am also going to derive an alternative expression that
I ended up using in mize. Finally, I will write some R code to
demonstrate their equivalence. This will turn out to be trickier than I
anticipated.

I originally tried to stick with the notation used by Sutskever but in
retrospect that was a mistake. It’s too hard to keep everything straight
between the different papers. I will borrow some symbols and try and be
very clear about what they mean and what each expression applies to.
Translating between this document and the original papers may require
mentally transposing some symbols and shifting iteration indices from
\\t\\ to \\t+1\\ or \\t-1\\.

Apart from the Sutskever and Bengio papers linked to above, it’s worth
checking out the [appendix to the Sutskever paper
(PDF)](https://proceedings.mlr.press/v28/sutskever13-supp.pdf) or the
relevant part (chapter 7) of [Sutskever’s
thesis](https://hdl.handle.net/1807/36012). For an introduction to NAG
itself, try the first part of this paper by [O’Donoghue and
Candès](https://arxiv.org/abs/1204.3982) (the rest of it is also good
but not germane to Nesterov momentum expressions).

## Definitions

The goal is to optimize some parameters, \\\theta\\. This is going to
involve gradient descent, so we will be evaluating the gradient of an
objective function of those parameters, \\\nabla f\left(\theta\right)\\,
and moving a certain distance in the direction of the negative of the
gradient, the distance being related to the learning rate,
\\\varepsilon\\. There will also be a momentum term involved, with
momentum coefficient \\\mu\\. The value of the parameters, learning rate
and momentum at iteration \\t\\ will be indicated by a subscript,
e.g. \\\theta_t\\.

### The Update Vector

A definition which holds across all the methods discussed here is that
the parameters at iteration \\t+1\\ are related to the parameters at
iteration \\t\\ by an update that involves the addition of a velocity
vector, \\v\\:

\\\theta\_{t+1} = \theta_t + v\_{t+1} \implies v\_{t+1} =
\theta\_{t+1} - \theta_t\\

This is pretty obvious, but the expressions we will be dealing with here
use a recursive definition of \\v_t\\, i.e. to usefully define these
updates we need to express \\v_t\\ in terms of \\v\_{t-1}\\. The exact
expression depends on whether we are using classical momentum or NAG. In
both cases, it would make sense to label a particular vector as \\v_t\\,
but due to the similarity in how they work, the temptation to substitute
the classical momentum recursive definition of \\v_t\\ into an
expression for NAG is almost overwhelming, and also always wrong.

To confuse matters further, some papers (see the Bengio paper for
example) directly use \\v_t\\ in expressions that update parameters, but
where the \\v_t\\ doesn’t actually apply to the parameter being updated.

Finally, I’ll try to be strict about names. When I’m referring to the
original Nesterov Accelerated Gradient, I’ll call that NAG. When I’m
referring to the versions that recast it as a momentum scheme, I’ll call
that Nesterov momentum, and refer to either the Sutskever or Bengio
formulation where necessary. The traditional momentum scheme is referred
to as “standard” or “regular” momentum by Bengio, and “classical” by
Sutskever. Other authors call it “Polyak” or “heavy ball” momentum. I’ll
refer to it as “classical” momentum.

## Classical Momentum

The classical momentum velocity vector is defined by Bengio as:

\\ v_t = \mu\_{t-1}v\_{t-1} - \varepsilon\_{t-1}\nabla
f\left(\theta\_{t-1}\right) \\

Sutskever gives the same definition but without any \\t\\ subscript on
the velocity vector or the momentum coefficient.

We can write out the full update for classical momentum as:

\\\theta\_{t+1} = \theta_t + \mu_t v_t - \varepsilon_t \nabla
f\left(\theta_t\right)\\

with velocity:

\\v\_{t+1} = \mu_t v_t - \varepsilon_t \nabla f\left(\theta_t \right)\\

Later we will want to use the recursive definition of momentum, but it
won’t be the full update, so we will not be able to define it as
\\v_t\\. In those circumstances, I will call the update \\m_t\\ to avoid
confusion:

\\ m\_{t+1} = \mu_t m_t - \varepsilon_t \nabla f\left(\theta_t \right)
\\

Remember: \\v\_{t+1} = m\_{t+1}\\ *only* in the case of classical
momentum.

## NAG

The Nesterov Accelerated Gradient method consists of a gradient descent
step, followed by something that looks a lot like a momentum term, but
isn’t exactly the same as that found in classical momentum. I’ll call it
a “momentum stage” here. If you look at the Sutskever paper, It’s
important to note that the parameters being minimized by NAG are given
the symbol \\y\\, not \\\theta\\. You’ll see that \\\theta\\ is the
symbol for the parameters after they’ve been updated by the gradient
descent stage, but before the momentum stage.

I’m going to stick with \\\theta_t\\ as the parameters being updated.
But we do need to define a symbol that accounts for the state of the
parameters after the gradient descent stage. I’ll go with
\\\phi\_{t+1}\\. Whether you consider the intermediate step to be part
of iterate \\t\\ or \\t+1\\ doesn’t really matter, but different papers
use different conventions, so you will need to keep track of it all or
risk getting very confused. In this case I am going with \\t+1\\ so it’s
consistent with any new quantities calculated in an iteration (like
\\v\_{t+1}\\).

Here’s the gradient descent stage:

\\ \phi\_{t+1} = \theta_t - \varepsilon_t \nabla f\left(\theta_t \right)
\\

followed by the momentum-like stage:

\\ \theta\_{t+1} = \phi\_{t+1} + \mu_t\left(\phi\_{t+1} - \phi\_{t}
\right) \\

That concludes one iteration of NAG. The hard part is choosing the
learning rate and momentum to get the convergence guarantees that make
the method attractive. Here I’ll stick to comparing the update
expressions, using fixed values in the examples later on to keep things
simple.

## Sutskever Nesterov Momentum

The key idea behind the Sutskever momentum derivation is to shift the
perspective about which of the parameters we want as the result of the
iteration, from \\\theta\\ to \\\phi\\. Rather than having the
optimization iterations proceed as “gradient descent, momentum (end of
iteration 1), gradient descent, momentum (end of iteration 2), gradient
descent etc.” move the boundary of where the iterations end by a
half-iteration to get “momentum, gradient descent (end of iteration 1),
momentum, gradient descent (end of iteration 2) etc.”. This leaves a
phantom gradient descent step that used to be first stage of the first
iteration, now floating in the nether regions of iteration zero, but you
can just pretend that the starting position is the result of gradient
descent from some other arbitrary starting position.

Perhaps a visualization will help. Here are two full iterations of NAG
in its standard form, separated by boxes so you can see where one
iteration ends and the other begins.

\\ \begin{aligned} \boxed{\begin{aligned} \phi\_{t+1} &= \theta_t -
\varepsilon_t \nabla f\left(\theta_t \right) \\ \theta\_{t+1} &=
\phi\_{t+1} + \mu_t\left(\phi\_{t+1} - \phi\_{t} \right) \end{aligned}}
\\\[0.75em\] \boxed{\begin{aligned} \phi\_{t+2} &= \theta\_{t+1} -
\varepsilon\_{t+1} \nabla f\left(\theta\_{t+1} \right) \\ \theta\_{t+2}
&= \phi\_{t+2} + \mu\_{t+1}\left(\phi\_{t+2} - \phi\_{t+1} \right)
\end{aligned}} \end{aligned} \\

The Sutskever derivation says shift the boundaries of the iteration so
that it starts one stage later. This is the same four equations as
above, but with a box showing the new definition of the iteration:

\\ \begin{gathered} \phi\_{t+1} = \theta_t - \varepsilon_t \nabla
f\left(\theta_t \right) \\ \boxed{\begin{aligned} \theta\_{t+1} &=
\phi\_{t+1} + \mu_t\left(\phi\_{t+1} - \phi\_{t} \right) \\ \phi\_{t+2}
&= \theta\_{t+1} - \varepsilon\_{t+1} \nabla f\left(\theta\_{t+1}
\right) \end{aligned}} \\ \theta\_{t+2} = \phi\_{t+2} +
\mu\_{t+1}\left(\phi\_{t+2} - \phi\_{t+1} \right) \end{gathered} \\

Now that the definition of where we start and finish each iteration has
changed, let’s shift all the indices in the box back one step. That
includes the momentum coefficient, which becomes \\\mu\_{t-1}\\:

\\ \theta_t = \phi_t + \mu\_{t-1}\left(\phi_t - \phi\_{t-1} \right) \\
\phi\_{t+1} = \theta_t - \varepsilon_t \nabla f\left(\theta_t \right) \\

Then because our update is now in terms of \\\phi\\ rather than
\\\theta\\, we can label \\\phi_t - \phi\_{t-1}\\ as \\v_t\\. The
lookahead point is then:

\\ \theta_t = \phi_t + \mu\_{t-1} v_t \\

And we can now substitute that definition of \\\theta_t\\ into the
gradient stage, which lets us write this half-shifted NAG iteration in
one line:

\\ \phi\_{t+1} = \phi_t + \mu\_{t-1} v_t - \varepsilon_t \nabla
f\left(\phi_t + \mu\_{t-1} v_t \right) \\

If you are comparing this with a version that uses \\\mu_t\\, remember
to shift the momentum schedule along with the iteration numbers. With
constant momentum, that’s one less thing to worry about.

There you have it: this looks just like the classical momentum update,
except that the gradient is calculated at \\\phi_t + \mu\_{t-1}v_t\\,
after the momentum update. Classical momentum would evaluate it at
\\\phi_t\\. Hence, one can do NAG by simply reversing the order in which
the update is usually carried out: do the momentum stage first, update
the parameters, and then do the gradient descent part. Just like
Sutskever said you could.

## Bengio Nesterov Momentum

I found working through the Bengio version of Nesterov momentum a bit
brain-breaking. It starts from the Sutskever definition and then defines
a new variable, \\\Theta\\, which represents the parameters after the
momentum update. But that just moves us one half-step forward in the
iteration, so we are back to doing things in the original NAG way:
gradient descent then momentum-style update. So their paper’s \\\Theta\\
is the same as our \\\theta\\. Its explanation of “committing to the
‘peekedahead’ parameters” and then “backtracking by the same amount
before each update”, also never landed with me (this is a reflection on
my mental capacity, not the paper). Finally, it expresses the update in
terms of \\v\\, but that’s still the the update in terms of the
difference in \\\phi\\, not \\\theta\\. Unlike this document it also
uses \\t-1\\ and \\t\\ to index the parameters that are being updated,
so you will have to remember to add or subtract 1 from the subscripts
when looking between here and the paper, but that’s the least of your
problems honestly.

It’s all too much. I can’t get their explanation to stick in my brain,
even though I can follow the individual steps and get the right
expression. So I will show you how to arrange the standard NAG
expressions into their form, rather than starting from the Sutskever
result.

First, let’s remove the confusion around using \\v\\ as an update vector
which isn’t the change in the parameters that we are updating each
iteration. Instead of \\\phi\_{t+1} - \phi\_{t} = v\_{t+1}\\, define:

\\ \phi\_{t+1} - \phi\_{t} = b\_{t+1} \\

(the ‘b’ is for Bengio). Our goal now is to end up with an expression
that is defined using \\\theta_t\\ (no \\\phi_t\\ allowed) and \\b_t\\.

I’ll use \\s_t\\ as shorthand for the gradient step:

\\ s_t = -\varepsilon_t \nabla f\left(\theta_t\right). \\

The negative sign and the learning rate are both part of \\s_t\\.

Now we are going to go back to using the original NAG expression so to
help you scrub all the shifts in \\t\\ from your mind here is the
expression again:

\\ \phi\_{t+1} = \theta_t + s_t \\ \theta\_{t+1} = \phi\_{t+1} +
\mu_t\left(\phi\_{t+1} - \phi\_{t} \right) \\ Start by replacing all
uses of \\\phi\_{t+1}\\ in the second equation with the definition in
terms of \\\theta_t\\ given in the first equation:

\\ \begin{equation} \begin{split} \theta\_{t+1} & = \theta_t + s_t +
\mu_t \left( \theta_t + s_t - \phi\_{t} \right) \\ & = \left(1 + \mu_t
\right) \theta_t + \left(1 + \mu_t \right) s_t - \mu_t \phi_t
\end{split} \end{equation} \\

So far, so unpromising. Time to get \\b_t\\ involved. Remember how we
defined:

\\ \begin{equation} \begin{split} \theta\_{t+1} & = \phi\_{t+1} +
\mu_t\left(\phi\_{t+1} - \phi\_{t} \right) \\ & = \phi\_{t+1} + \mu_t
b\_{t+1} \end{split} \end{equation} \\

Let’s shift everything back one iteration:

\\ \theta\_{t} = \phi\_{t} + \mu\_{t-1} b\_{t} \\ \implies \phi\_{t} =
\theta\_{t} - \mu\_{t-1} b\_{t} \\

And now substitute for \\\phi_t\\ in the last term in our NAG update:

\\ \begin{equation} \begin{split} \theta\_{t+1} & = \left(1 + \mu_t
\right) \theta_t + \left(1 + \mu_t \right) s_t - \mu_t \phi_t \\ & =
\left(1 + \mu_t \right) \theta_t + \left(1 + \mu_t \right) s_t -
\mu\_{t} \theta_t + \mu\_{t-1} \mu_t b_t \end{split} \end{equation} \\

The final step is to gather together the \\\theta_t\\ terms and the
result is:

\\ \theta\_{t+1} = \theta_t + \mu\_{t-1} \mu_t b_t + \left(1 + \mu_t
\right) s_t \\

Ta da. There it is, just like in the paper, subject to the differences
in symbols and when the iteration ticks over from \\t\\ to \\t+1\\ and
so on.

The advantage of this expression for the Nesterov momentum is that it
doesn’t require calculating a gradient at a non-standard position, and
only requires a modification to the coefficients used to calculate the
velocity, which is probably an easier change to make to an existing
codebase which already uses classical momentum. However, because this
expression is recursive on \\b_t\\ rather than the update vector
\\v_t\\, you will need to remember to store a different vector than you
would with classical momentum. It’s not really clear to me whether that
is an easy or hard thing to do with most optimization code bases vs
having to make code changes to calculate the gradient at a different
location than the current parameter vector, as you have to do with the
Sutskever version.

You should also remember that Bengio gives us the parameters after the
momentum stage, while Sutskever gives us the parameters after gradient
descent. They’re half a step out of sync. To get matching results, we’ll
need to compare the same stage and get the starting vectors and momentum
coefficients lined up. This will occupy us for quite a while later on.

## An Alternative Expression for NAG

Let’s return to the two stages of NAG, keeping our original notation:
\\\phi\_{t+1}\\ is the position after gradient descent, and
\\\theta\_{t+1}\\ is the position after the momentum stage. This time
we’ll write the update in terms of the actual displacement of
\\\theta\\, rather than the Bengio buffer.

\\\phi\_{t+1} = \theta_t + s_t\\

\\\theta\_{t+1} = \phi\_{t+1} + \mu_t\left(\phi\_{t+1} - \phi_t
\right)\\

Now, let’s just write out the momentum stage in terms of \\\theta\\,
substituting \\\phi\\ wherever we find it:

\\ \theta\_{t+1} = \theta_t + s_t + \mu_t \left\[\theta_t + s_t -
\theta\_{t-1} - s\_{t-1}\right\] \\

Rearranging:

\\ \theta\_{t+1} = \theta_t + \mu_t \left\[\theta_t - \theta\_{t-1} -
s\_{t-1} + s_t\right\] + s_t \\

Finally, we can substitute in \\v_t\\ for the first two terms in the
square brackets, to give:

\\ \theta\_{t+1} = \theta_t + \mu_t \left\[v_t - s\_{t-1} +
s_t\right\] + s_t \\

with velocity:

\\ v\_{t+1} = \mu_t \left\[v_t - s\_{t-1} + s_t\right\] + s_t \\

This looks a lot like the classical momentum expression, but with the
velocity vector modified to first remove the contribution of the
gradient descent from the previous iteration, and replace it with the
gradient descent contribution from the *current* iteration. Gives an
interesting insight into the idea of the Nesterov momentum using a form
of “lookahead” with the gradient descent.

You could also choose to expand the velocity expression to make it look
a bit like the Bengio formulation:

\\ v\_{t+1} = \mu_t \left\[v_t - s\_{t-1}\right\] + \left(1 +
\mu_t\right) s_t \\

but as this version can’t be expressed as the classical momentum form
with different coefficients, the way the Bengio formulation can be, it
probably doesn’t gain you anything in terms of implementation, except
you can expand it and rearrange it further to give:

\\ v\_{t+1} = \mu_t v_t + s_t + \mu_t\left\[s_t - s\_{t-1}\right\] \\
which now resembles classical momentum with an extra correction term.
For a fixed learning rate, this is the form in equation (1.4) of [Shi
and
co-workers](https://link.springer.com/article/10.1007/s10107-021-01681-8),
who call the extra term a “gradient correction”. User ‘denis’ also uses
this expression in an [answer on the Cross Validated Stack
Exchange](https://stats.stackexchange.com/a/233430), calling it
“gradient momentum”. I’ll use “gradient correction” here.

*December 25 2025*: Here’s a possible insight into the meaning of the
gradient correction.

The correction is the change in the gradient step, scaled by \\\mu_t\\.
When the steps are similar in both direction and size, \\s_t - s\_{t-1}
\approx 0\\, the update looks a lot like classical momentum. But
consider a steep ravine where we are bouncing from one side to the
other. If the step roughly reverses each iteration, \\s_t \approx
-s\_{t-1}\\, the last term becomes approximately \\2\mu_t s_t\\. This
gives the current gradient step more say in where we go next. Remember
that \\s_t\\ includes the learning rate, so changing that will also
affect this term.

### What happens on a quadratic?

Let’s start with just one parameter, so our quadratic objective is an
ordinary function:

\\ f(\theta) = \frac{1}{2}a\theta^2 + c\theta. \\

Here \\a\\ and \\c\\ are constants. Taking the derivative gives:

\\ f'(\theta) = a\theta + c. \\

The gradient symbol \\\nabla f\\ we’ve been using is just this
derivative when there is only one parameter. So the difference between
the gradients at two successive positions is:

\\ \begin{aligned} f'(\theta_t) - f'(\theta\_{t-1}) &= (a\theta_t + c) -
(a\theta\_{t-1} + c) \\ &= a(\theta_t - \theta\_{t-1}) \\ &= a v_t.
\end{aligned} \\

The last line uses our definition of the velocity, \\v_t = \theta_t -
\theta\_{t-1}\\. With a fixed learning rate, \\s_t = -\varepsilon
f'(\theta_t)\\, so:

\\ s_t - s\_{t-1} = -\varepsilon\left\[f'(\theta_t) -
f'(\theta\_{t-1})\right\] = -\varepsilon a v_t. \\

Now keep the momentum coefficient \\\mu\\ fixed too, and substitute into
our NAG expression:

\\ \begin{aligned} v\_{t+1} &= \mu v_t + s_t + \mu(s_t - s\_{t-1}) \\ &=
\mu v_t - \varepsilon f'(\theta_t) - \mu\varepsilon a v_t \\ &= \mu(1 -
\varepsilon a)v_t - \varepsilon f'(\theta_t). \end{aligned} \\

This looks like classical momentum with a momentum coefficient of
\\\mu(1-\varepsilon a)\\. And \\a=f''(\theta)\\ is the curvature of our
quadratic. NAG keeps nearly all its momentum where the curvature is
close to zero, and reduces it as the positive curvature increases,
reaching zero at \\\varepsilon a=1\\. Beyond that, the momentum
contribution reverses direction. Negative curvature increases the
momentum instead.

This result generalizes to quadratic objectives with multiple
parameters; see section A.2 of the [supplement to Sutskever et
al](https://proceedings.mlr.press/v28/sutskever13-supp.pdf) for the
derivation.

We can also choose how much of this correction to apply, which brings us
to unified momentum.

## Unified Momentum

The earliest example of unified momentum I’ve found is in a 2016 paper
by [Yang, Lin and Li](https://arxiv.org/abs/1604.03257).

The idea is to introduce an extra parameter, \\\lambda_t\\, that lets us
choose how much of the gradient step correction to use:

\\ v\_{t+1} = \mu_t v_t + s_t + \lambda_t\mu_t(s_t-s\_{t-1}) \\ Setting
\\\lambda_t = 0\\ gives classical momentum and \\\lambda_t = 1\\ gives
NAG.

For our one-parameter quadratic, the same substitution we used above
gives:

\\ v\_{t+1} = \mu(1 - \lambda\varepsilon a)v_t - \varepsilon
f'(\theta_t). \\

So \\\lambda\\ controls how much the curvature changes the momentum
coefficient.

### Storing the NAG update

The gradient correction seems to require storing the previous gradient
step as well as the velocity. For NAG, where \\\lambda_t=1\\, we can
combine them into one saved vector:

\\ r_t = v_t - s\_{t-1}. \\

This has a meaning in terms of the positions we’ve already defined.
Since \\\phi_t=\theta\_{t-1}+s\_{t-1}\\, after a completed update we
have:

\\ \begin{aligned} r_t &= (\theta_t-\theta\_{t-1})-s\_{t-1} \\ &=
\theta_t-\phi_t \\ &= \mu\_{t-1}b_t. \end{aligned} \\

So the saved vector is the gap between the gradient descent and momentum
positions. Substituting \\v_t-s\_{t-1}=r_t\\ into the NAG update gives
\\v\_{t+1}=s_t+\mu_t(r_t+s_t)\\. Removing \\s_t\\ to prepare the saved
vector for the next iteration leaves us with:

\\ \begin{aligned} s_t &= -\varepsilon_t\nabla f(\theta_t), \\ r\_{t+1}
&= \mu_t(r_t+s_t), \\ \theta\_{t+1} &= \theta_t+s_t+r\_{t+1}.
\end{aligned} \\

Only \\r_t\\ needs to persist alongside the parameters; the actual
displacement is \\v\_{t+1}=s_t+r\_{t+1}\\. These equations also work
when the learning rate or momentum changes.

Start with \\r_0=0\\, corresponding to \\\theta_0=\phi_0\\. The first
displacement is then \\(1+\mu_0)s_0\\. Setting \\\mu_0=0\\ gives an
ordinary gradient step; using nonzero momentum gives a longer first
step. We’ll return to that choice when we unroll the updates.

## Unrolling NAG and Classical Momentum

As a way to understand the difference between classical momentum and
NAG, let’s write out the first few steps of the optimization,
substituting in the recursive definitions with the previous step, and
see what emerges.

I’ll assume that the momentum is constant, so we can write \\\mu\\
without a subscript.

Initial coordinates are at \\\theta_0\\.

### Classical Momentum

With constant momentum, the classical momentum update is:

\\ \theta\_{t+1} = \theta_t + s_t + \mu v_t \\

And the first four updates (including steepest descent on the first
iteration) are as follows:

\\ \begin{equation} \begin{split} \theta_1 & = \theta_0 + s_0 \\ \\
\theta_2 & = \theta_1 + s_1 + \mu v_1 \\ & = \theta_1 + s_1 + \mu s_0 \\
\\ \theta_3 & = \theta_2 + s_2 + \mu v_2 \\ & = \theta_2 + s_2 + \mu
\left( s_1 + \mu s_0 \right) \\ & = \theta_2 + s_2 + \mu s_1 + \mu^2 s_0
\\ \\ \theta_4 & = \theta_3 + s_3 + \mu v_3 \\ & = \theta_3 + s_3 + \mu
\left( s_2 + \mu s_1 + \mu^2 s_0 \right)\\ & = \theta_3 + s_3 + \mu
s_2 + \mu^2 s_1 + \mu^3 s_0 \\ \end{split} \end{equation} \\

### NAG

We’ll use this straight-forward expression for NAG:

\\ \theta\_{t+1} = \theta_t + s_t + \mu \left(\theta_t + s_t -
\theta\_{t-1} - s\_{t-1} \right) \\

with just this bit of re-arranging so terms are all grouped together and
the \\v_t\\ is at the end, which makes it easier to follow the
substitution of \\v_t\\:

\\ \theta\_{t+1} = \theta_t + \left(1 + \mu \right) s_t - \mu s\_{t-1} +
\mu v_t \\

#### The first step

*December 14 2021*: This section describing your choices for the first
step is new.

Like classical momentum, we have no previous directions to draw on for
the first step. So let’s start both sets of parameters at the same
place, \\\theta_0=\phi_0\\, and see what happens. The gradient stage
gives us:

\\ \phi_1 = \theta_0 + s_0, \\

and then the momentum stage gives us:

\\ \theta_1 = \phi_1 + \mu_0(\phi_1-\phi_0) = \theta_0 + (1+\mu_0)s_0.
\\

There’s the longer first step: with nonzero momentum, we move
\\(1+\mu_0)\\ times as far along the negative gradient. The momentum
stage already has something to work with, because \\\phi_1-\phi_0\\ is
the gradient step we just took.

If we set \\\mu_0=0\\, the two stages finish at the same point and we
get the same first step as classical momentum. This is what the schedule
in the [Sutskever
supplement](https://proceedings.mlr.press/v28/sutskever13-supp.pdf)
does. We’ll look at both choices, because this apparently minor detail
will turn up again when we try to get the code to agree.

First, here’s the “long” version, using the same nonzero \\\mu\\ from
the start and taking our results after the momentum stage:

\\ \begin{equation} \begin{split} \theta_1 & = \theta_0 + \left( 1 + \mu
\right) s_0 \\ \\ \theta_2 & = \theta_1 + \left( 1 + \mu \right) s_1 -
\mu s_0 + \mu v_1 \\ & = \theta_1 + \left(1 + \mu \right) s_1 - \mu
s_0 + \mu \left( 1 + \mu \right) s_0 \\ & = \theta_1 + \left(1 + \mu
\right) s_1 + \mu^2 s_0 \\ \\ \theta_3 & = \theta_2 + \left( 1 + \mu
\right) s_2 - \mu s_1 + \mu v_2 \\ & = \theta_2 + \left(1 + \mu \right)
s_2 - \mu s_1 + \mu \left(1 + \mu \right) s_1 + \mu^3 s_0\\ & =
\theta_2 + \left(1 + \mu \right) s_2 + \mu^2 s_1 + \mu^3 s_0 \\ \\
\theta_4 & = \theta_3 + \left(1 + \mu \right) s_3 - \mu s_2 + \mu v_3 \\
& = \theta_3 + \left(1 + \mu \right) s_3 - \mu s_2 + \mu \left\[
\left(1 + \mu \right) s_2 + \mu^2 s_1 + \mu^3 s_0 \right\] \\ & =
\theta_3 + \left(1 + \mu \right) s_3 - \mu s_2 + \mu \left(1 + \mu
\right) s_2 + \mu^3 s_1 + \mu^4 s_0 \\ & = \theta_3 + \left(1 + \mu
\right) s_3 + \mu^2 s_2 + \mu^3 s_1 + \mu^4 s_0 \end{split}
\end{equation} \\

You can see a pattern emerging here at this point, I hope. These are
quite similar in form compared to the classical momentum.

### NAG with a “short” first step

What happens if we set \\\mu_0=0\\ and then use the fixed value \\\mu\\
from the second iteration onwards? We start with the same gradient step
as classical momentum. Here’s how that plays out:

\\ \begin{equation} \begin{split} \theta\_{1} & = \theta_0 + s_0 \\ \\
\theta_2 & = \theta_1 + \left( 1 + \mu \right) s_1 - \mu s_0 + \mu v_1
\\ & = \theta_1 + \left( 1 + \mu \right) s_1 - \mu s_0 + \mu s_0 \\ & =
\theta_1 + \left( 1 + \mu \right) s_1 \\ \\ \theta_3 & = \theta_2 +
\left( 1 + \mu \right) s_2 - \mu s_1 + \mu v_2 \\ & = \theta_2 +
\left(1 + \mu \right) s_2 - \mu s_1 + \mu \left(1 + \mu \right) s_1 \\ &
= \theta_2 + \left(1 + \mu \right) s_2 + \mu^2 s_1 \\ \\ \theta_4 & =
\theta_3 + \left( 1 + \mu \right) s_3 - \mu s_2 + \mu v_3 \\ & =
\theta_3 + \left( 1 + \mu \right) s_3 - \mu s_2 + \mu \left\[ \left( 1 +
\mu \right) s_2 + \mu^2 s_1 \right\] \\ & = \theta_3 + \left( 1 + \mu
\right) s_3 - \mu s_2 + \mu \left( 1 + \mu \right) s_2 + \mu^3 s_1 \\ &
= \theta_3 + \left( 1 + \mu \right) s_3 + \mu^2 s_2 + \mu^3 s_1
\end{split} \end{equation} \\

Overall the results aren’t that different: nearly all the terms have the
same coefficients. But the \\s_0\\ terms cancel out of the second
update, and once they’re gone, they’re gone for good. The first gradient
still got us to \\\theta_1\\, of course; we just stop adding a
contribution from it in later updates. And the second update is
\\(1+\mu)s_1\\, so the “long” gradient step has turned up again, one
iteration later.

#### Relative Weights

Here’s a table showing how, on the fourth update from the initial point,
the methods weight the current and previous gradient steps:

|         Momentum Type          |   \\s_3\\   |  \\s_2\\  |  \\s_1\\  |  \\s_0\\  |
|:------------------------------:|:-----------:|:---------:|:---------:|:---------:|
|           Classical            |      1      |  \\\mu\\  | \\\mu^2\\ | \\\mu^3\\ |
| NAG, fixed \\\mu\\ from step 0 | 1 + \\\mu\\ | \\\mu^2\\ | \\\mu^3\\ | \\\mu^4\\ |
|    NAG, \\\mu_0=0\\ startup    | 1 + \\\mu\\ | \\\mu^2\\ | \\\mu^3\\ |     0     |

With fixed \\\mu\\, NAG gives the current gradient step more weight and
multiplies the coefficients of all the older steps by an extra factor of
\\\mu\\. After that, the weights decay at the same rate for both
methods: each older step gets another factor of \\\mu\\.

There’s a nice bit of bookkeeping here too. For \\0\leq\mu\<1\\, if we
keep adding older terms, the NAG weights sum to:

\\ 1 + \mu + \frac{\mu^2}{1-\mu} = \frac{1}{1-\mu}. \\

That’s the same sum as classical momentum. NAG has shifted some weight
from the older steps to the current one, while keeping the total the
same.

Dividing the newest coefficient by this sum gives a long-history weight
of \\1-\mu\\ for classical momentum and \\(1+\mu)(1-\mu)=1-\mu^2\\ for
NAG. At \\\mu=0.9\\, that’s 10% versus 19%. The four-step examples below
have much less history to share the weight with, so their percentages
are larger.

The regularity of the difference between classical momentum and NAG
weights suggest that there should be a way to express the NAG update in
terms of a classical momentum update. See the section on ‘Dozat Nesterov
momentum’ below for more on that.

Let’s write some R code to generate a table to compare how the relative
weights turn out between CM and NAG at different values of \\\mu\\.
We’ll normalize the coefficients so they sum to 1.

``` r

cm_weights <- function(mu) { c(1, mu, mu * mu, mu * mu * mu) }

nag_zero_start_weights <- function(mu) { c(1 + mu, mu * mu, mu * mu * mu, 0) }

nag_weights <- function(mu) { c(1 + mu, mu * mu, mu * mu * mu, mu * mu * mu * mu) }


mumat <- function(wfun, mus = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)) {
  rel_weights <- matrix(nrow = length(mus), ncol = 5)
  for (i in 1:length(mus)) {
    mu <- mus[i]
    weights <- wfun(mu)
    rel_weights[i, ] <- c(mu, weights / sum(weights))
  }
  colnames(rel_weights) <- c("mu", "s3", "s2", "s1", "s0")

  rel_weights
}
```

#### Classical momentum weights

``` r

knitr::kable(mumat(cm_weights), digits = 4)
```

|   mu |     s3 |     s2 |     s1 |     s0 |
|-----:|-------:|-------:|-------:|-------:|
| 0.10 | 0.9001 | 0.0900 | 0.0090 | 0.0009 |
| 0.25 | 0.7529 | 0.1882 | 0.0471 | 0.0118 |
| 0.50 | 0.5333 | 0.2667 | 0.1333 | 0.0667 |
| 0.75 | 0.3657 | 0.2743 | 0.2057 | 0.1543 |
| 0.90 | 0.2908 | 0.2617 | 0.2355 | 0.2120 |
| 0.99 | 0.2538 | 0.2512 | 0.2487 | 0.2462 |

#### NAG weights

``` r

knitr::kable(mumat(nag_weights), digits = 4)
```

|   mu |     s3 |     s2 |     s1 |     s0 |
|-----:|-------:|-------:|-------:|-------:|
| 0.10 | 0.9900 | 0.0090 | 0.0009 | 0.0001 |
| 0.25 | 0.9384 | 0.0469 | 0.0117 | 0.0029 |
| 0.50 | 0.7742 | 0.1290 | 0.0645 | 0.0323 |
| 0.75 | 0.5736 | 0.1844 | 0.1383 | 0.1037 |
| 0.90 | 0.4640 | 0.1978 | 0.1780 | 0.1602 |
| 0.99 | 0.4060 | 0.2000 | 0.1980 | 0.1960 |

#### NAG weights with \\\mu_0=0\\ startup

``` r

knitr::kable(mumat(nag_zero_start_weights), digits = 4)
```

|   mu |     s3 |     s2 |     s1 |  s0 |
|-----:|-------:|-------:|-------:|----:|
| 0.10 | 0.9901 | 0.0090 | 0.0009 |   0 |
| 0.25 | 0.9412 | 0.0471 | 0.0118 |   0 |
| 0.50 | 0.8000 | 0.1333 | 0.0667 |   0 |
| 0.75 | 0.6400 | 0.2057 | 0.1543 |   0 |
| 0.90 | 0.5525 | 0.2355 | 0.2120 |   0 |
| 0.99 | 0.5050 | 0.2487 | 0.2462 |   0 |

For the four steps shown here, at high momentum CM puts around 25–30% of
the weight on the current step, compared with 41–46% for NAG and 50–55%
for NAG with the short first step. That’s quite a shift in emphasis,
even after just four iterations.

### Unrolling with a momentum schedule

I kept \\\mu\\ constant in the unrolling expressions so far to keep
things neat and to highlight some similarities between CM and NAG. But
what if we want \\\mu\\ to vary with iteration? In our derivations of
Nesterov momentum schemes we have explicitly accounted for \\\mu_t\\ so
let’s take a look (and brace ourselves for the expressions to lose a lot
of readability):

#### Classical Momentum with a Schedule

\\ \begin{equation} \begin{split} \theta_1 & = \theta_0 + s_0 \\ \\
\theta_2 & = \theta_1 + s_1 + \mu_1 v_1 \\ & = \theta_1 + s_1 + \mu_1
s_0 \\ \\ \theta_3 & = \theta_2 + s_2 + \mu_2 v_2 \\ & = \theta_2 +
s_2 + \mu_2 \left( s_1 + \mu_1 s_0 \right) \\ & = \theta_2 + s_2 + \mu_2
s_1 + \mu_2 \mu_1 s_0 \\ \\ \theta_4 & = \theta_3 + s_3 + \mu_3 v_3 \\ &
= \theta_3 + s_3 + \mu_3 \left( s_2 + \mu_2 s_1 + \mu_2 \mu_1 s_0
\right)\\ & = \theta_3 + s_3 + \mu_3 s_2 + \mu_3 \mu_2 s_1 + \mu_3 \mu_2
\mu_1 s_0 \\ \end{split} \end{equation} \\ Already this is looking
pretty gross.

#### NAG with a Momentum Schedule

\\ \begin{equation} \begin{split} \theta_1 & = \theta_0 + \left( 1 +
\mu_0 \right) s_0 \\ \\ \theta_2 & = \theta_1 + \left( 1 + \mu_1 \right)
s_1 - \mu_1 s_0 + \mu_1 v_1 \\ & = \theta_1 + \left(1 + \mu_1 \right)
s_1 - \mu_1 s_0 + \mu_1 \left( 1 + \mu_0 \right) s_0 \\ & = \theta_1 +
\left(1 + \mu_1 \right) s_1 + \mu_1 \mu_0 s_0 \\ \\ \theta_3 & =
\theta_2 + \left( 1 + \mu_2 \right) s_2 - \mu_2 s_1 + \mu_2 v_2 \\ & =
\theta_2 + \left(1 + \mu_2 \right) s_2 - \mu_2 s_1 + \mu_2 \left(1 +
\mu_1 \right) s_1 + \mu_2 \mu_1 \mu_0 s_0\\ & = \theta_2 + \left(1 +
\mu_2 \right) s_2 + \mu_2 \mu_1 s_1 + \mu_2 \mu_1 \mu_0 s_0 \\ \\
\theta_4 & = \theta_3 + \left(1 + \mu_3 \right) s_3 - \mu_3 s_2 + \mu_3
v_3 \\ & = \theta_3 + \left(1 + \mu_3 \right) s_3 - \mu_3 s_2 + \mu_3
\left\[ \left(1 + \mu_2 \right) s_2 + \mu_2 \mu_1 s_1 + \mu_2 \mu_1
\mu_0 s_0 \right\] \\ & = \theta_3 + \left(1 + \mu_3 \right) s_3 - \mu_3
s_2 + \mu_3 \left(1 + \mu_2 \right) s_2 + \mu_3 \mu_2 \mu_1 s_1 + \mu_3
\mu_2 \mu_1 \mu_0 s_0 \\ & = \theta_3 + \left(1 + \mu_3 \right) s_3 +
\mu_3 \mu_2 s_2 + \mu_3 \mu_2 \mu_1 s_1 + \mu_3 \mu_2 \mu_1 \mu_0 s_0
\end{split} \end{equation} \\

I can see why I was in no hurry to write all that out in full until now.
With a \\\mu_0=0\\ startup, the \\s_0\\ term vanishes after the first
returned point, just as in the constant-\\\mu\\ calculation above.

Here’s a comparison of the weights applied to the different gradient
steps at \\\theta_4\\:

| Momentum Type | \\s_3\\ | \\s_2\\ | \\s_1\\ | \\s_0\\ |
|:--:|:---|:---|:---|:---|
| Classical | 1 | \\\mu_3\\ | \\\mu_3 \mu_2\\ | \\\mu_3 \mu_2 \mu_1\\ |
| NAG | 1 + \\\mu_3\\ | \\\mu_3 \mu_2\\ | \\\mu_3 \mu_2 \mu_1\\ | \\\mu_3 \mu_2 \mu_1 \mu_0\\ |

Now the momentum is not constant at each iteration, we see a more
complex relationship between the CM and NAG weights. There definitely
doesn’t seem to be a simple re-weighting that could map from CM to NAG
here: each gradient needs re-weighting by a different \\\mu_t\\.

## Dozat Nesterov Momentum

A paper by [Dozat integrating Nesterov
momentum](https://openreview.net/forum?id=OM0jvwB8jIp57ZJjtNEZ) into the
[Adam](https://arxiv.org/abs/1412.6980) gives a slightly different
expression, that is quite closely related to the Bengio form:

\\ \begin{equation} \begin{split} \theta\_{t+1} & = \phi\_{t+1} +
\mu\_{t}\left(\phi\_{t+1} - \phi\_{t} \right) \\ & = \phi\_{t+1} +
\mu\_{t} b\_{t+1} \\ & = \theta\_{t} + s\_{t} + \mu\_{t} b\_{t+1}
\end{split} \end{equation} \\ Now let us express \\b\_{t+1}\\
recursively in terms of the previous \\b\_{t}\\:

\\ \begin{equation} \begin{split} b\_{t+1} & = \phi\_{t+1} - \phi\_{t}
\\ & = \left(\theta\_{t} + s\_{t} \right) - \phi\_{t} \\ & = \phi\_{t} +
\mu\_{t-1} \left(\phi\_{t} - \phi\_{t-1} \right) + s\_{t} - \phi\_{t} \\
& = \mu\_{t-1} b\_{t} + s\_{t} \end{split} \end{equation} \\ This looks
just like the classical momentum update, with \\b_t\\ instead of \\m_t\\
or \\v_t\\ and the momentum coefficient shifted back one iteration.
Let’s put the two forms next to each other.

Classical momentum looks like:

\\ \begin{equation} \begin{split} m\_{t+1} = \mu_t m\_{t} + s\_{t} \\
\theta\_{t+1} = \theta\_{t} + m\_{t+1} \end{split} \end{equation} \\ and
NAG looks like:

\\ \begin{equation} \begin{split} b\_{t+1} = \mu\_{t-1} b\_{t} + s\_{t}
\\ \theta\_{t+1} = \theta\_{t} + s\_{t} + \mu\_{t} b\_{t+1} \end{split}
\end{equation} \\

So a very similar recursive form for the momentum buffer for both CM and
NAG, but using momentum from a different time step, and then a different
velocity vector for the update. The Dozat paper shifts the momentum
indices by one time step, so remember to shift the schedule too when
translating between the two. Constant momentum saves us some trouble
here again.

### Relation to QHM

This way of combining a momentum buffer with the current gradient step
also turns up in [quasi-hyperbolic momentum
(QHM)](https://arxiv.org/abs/1810.06801), introduced by Ma and Yarats in
2018. I wrote some extra [notes on
QHM](https://jlmelville.github.io/mize/articles/qhm.html) for my own
benefit, but we can make the connection here with the expressions we
already have.

With fixed momentum, the buffer \\b_t\\ above follows the classical
momentum recurrence:

\\ b\_{t+1} = \mu b_t + s_t. \\

Now recall [unified momentum](#unified-momentum). Keeping \\\mu\\,
\\\varepsilon\\ and \\\lambda\\ fixed, with \\0\leq\mu\<1\\, its update
is:

\\ v\_{t+1} = \mu v_t + s_t + \lambda\mu(s_t - s\_{t-1}) \\

Starting with \\v_0=b_0=s\_{-1}=0\\, a bit of substitution gives the
following expression, which you can check by induction:

\\ \boxed{ v\_{t+1} = \lambda s_t + \left\[1 - \lambda(1 -
\mu)\right\]b\_{t+1}. } \\

Setting \\\lambda=1\\ gives \\v\_{t+1}=s_t+\mu b\_{t+1}\\, exactly the
Dozat update above. Setting \\\lambda=0\\ gives \\v\_{t+1}=b\_{t+1}\\,
the classical momentum update.

QHM normalizes its momentum buffer by a factor of \\1-\mu\\, which means
we also need to rescale the learning rate. Translating into the QHM
paper’s symbols gives:

\\ \beta\_{\mathrm{QHM}} = \mu, \qquad \nu\_{\mathrm{QHM}} = 1 -
\lambda(1 - \mu), \qquad \alpha\_{\mathrm{QHM}} = \frac{\varepsilon}{1 -
\mu}. \\

With \\\lambda=0\\, we get \\\nu\_{\mathrm{QHM}}=1\\, which is classical
momentum. With \\\lambda=1\\, we get \\\nu\_{\mathrm{QHM}}=\mu\\, which
is NAG. So letting \\\lambda\\ run from 0 to 1 takes us along the part
of QHM between classical momentum and NAG.

### Where the learning rate lives

There’s another implementation detail that changes the weighting of
history: does an old gradient keep the learning rate it had when we
calculated it, or does today’s learning rate rescale its contribution?
To isolate that question, keep \\\mu\\ fixed and write \\g_t=\nabla
f(\theta_t)\\.

Our buffer stores gradient steps, including their learning rates:

\\ \begin{aligned} b\_{t+1} &= \mu b_t-\varepsilon_t g_t, \\
\theta\_{t+1} &= \theta_t-\varepsilon_t g_t+\mu b\_{t+1}. \end{aligned}
\\

Alternatively, we can store the gradients in a buffer \\B_t\\, and apply
the learning rate when we update the parameters:

\\ \begin{aligned} B\_{t+1} &= \mu B_t+g_t, \\ \theta\_{t+1} &=
\theta_t-\varepsilon_t(g_t+\mu B\_{t+1}). \end{aligned} \\

With a constant learning rate and both buffers starting at zero,
\\b_t=-\varepsilon B_t\\ makes the updates equivalent. With a changing
learning rate, the first version preserves the old learning rates inside
its accumulated steps, while the second applies the current rate to the
whole gradient history. They can therefore give different results under
the same learning-rate schedule.

These are both used in practice:
[Keras](https://keras.io/api/optimizers/sgd/) documents the first
convention, while
[PyTorch](https://docs.pytorch.org/docs/stable/generated/torch.optim.SGD.html)
documents the second and notes the difference for Nesterov momentum too.
We’re comparing just the momentum updates here, with no weight decay or
dampening.

## NAG in practice

So that’s how it’s all supposed to work in principle. Here I’ll
demonstrate the equivalence of these methods, by implementing them all
in simple R code.

The biggest simplification I’ll make is that I’ll assume a constant
learning rate and a constant momentum coefficient.

### Classical Momentum

``` r

#' Optimization by Classical Momentum
#'
#' @param par Starting point of vector of parameters to optimize.
#' @param fn Objective function to optimize. Takes vector with length of
#' \code{par} and returns a scalar.
#' @param gr Gradient of the objective function \code{fn}. Takes vector with
#' length of \code{par} and returns a vector of the same length.
#' @param lr Learning rate.
#' @param mu Momentum coefficient.
#' @param max_iter Maximum number of iterations to optimize for. First iteration
#' is always steepest descent.
#' @return list with components: \code{par} final set of parameters; \code{f}
#' value of \code{fn} evaluated at the returned set of parameters; \code{fs}
#' vector of function evaluated after each iteration.
cm <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)

  v <- rep(0, length(par))
  for (i in 1:max_iter) {
    g <- gr(par)
    v <- mu * v - lr * g
    par <- par + v

    # store results
    f <- fn(par)
    fs[i] <- f
  }

  list(par = par, f = f, fs = fs)
}
```

This is a reference implementation of classical momentum. Nearly all the
parameters and the return values are the same for the other functions
(except where noted), so they’re documented here once.

Onto the various NAG implementations. The `nag` routine below sets the
momentum coefficient to zero on the first iteration, giving us the
“short” first step. The Bengio, Dozat and momentum-style versions use
the nonzero value straight away. Sutskever gives us its result after the
gradient stage. So we already have a few differences to keep an eye on.
Classical momentum gets its initial gradient step simply by starting the
velocity vector at zero.

### NAG

``` r

# Optimization by Nesterov Accelerated Gradient
#
# Return list also contains gd_fs: function evaluated after gradient descent
# stage of each iteration; all: function evaluated after gradient descent
# stage and momentum stage, in order.
nag <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)
  gd_fs <- rep(0, max_iter)
  all <- rep(0, max_iter * 2)

  x_old <- rep(0, length(par))
  for (i in 1:max_iter) {
    # gradient descent stage
    g <- gr(par)
    x <- par - (lr * g)

    # store gradient descent values
    f <- fn(x)
    gd_fs[i] <- f
    all[i * 2 - 1] <- f

    # momentum stage and update
    par <- x + ifelse(i == 1, 0, mu) * (x - x_old)
    x_old <- x

    # store momentum values
    f <- fn(par)
    fs[i] <- f
    all[i * 2] <- f
  }

  list(par = par, f = f, fs = fs, gd_fs = gd_fs, all = all)
}
```

In this routine, rather than store the velocity vector, we store the
previous gradient descent result, `x_old`. In order to ensure the first
iteration is gradient descent only, we also need to manually set the
momentum coefficient `mu` to zero on the first iteration, which is what
that `ifelse` expression does.

Also, there’s some extra code to calculate and store the function values
after the gradient descent stage. These aren’t needed for the
optimization to work, I just want to keep track of the values to
demonstrate the different methods are in fact equivalent.

### Sutskever Formulation

``` r

# Optimization by Sutskever Nesterov Momentum
#
# Return list also contains mu_fs: function evaluated after momentum
# stage of each iteration; all: function evaluated after gradient descent
# stage and momentum stage, in order.
snag <- function(par, fn, gr, lr, mu, max_iter = 10) {
  v <- rep(0, length(par))

  fs <- rep(0, max_iter)
  mu_fs <- rep(0, max_iter)
  all <- rep(0, max_iter * 2)

  for (i in 1:max_iter) {
    # momentum stage and update parameters
    mu_step <- mu * v
    par <- par + mu_step

    # store momentum results
    f <- fn(par)
    mu_fs[i] <- f
    all[i * 2 - 1] <- f

    # gradient descent stage
    g <- gr(par)
    gd_step <- -lr * g

    # update and store velocity for next step
    par <- par + gd_step
    v <- mu_step + gd_step

    # store gradient descent results
    f <- fn(par)
    fs[i] <- f
    all[i * 2] <- f
  }

  list(par = par, f = f, fs = fs, mu_fs = mu_fs, all = all)
}
```

This one is pretty straight-forward.

### Bengio Formulation

``` r

# Optimization by Bengio Nesterov Momentum
bnag <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)

  v <- rep(0, length(par))
  for (i in 1:max_iter) {
    g <- gr(par)
    # here mu_{t-1} * mu_t equals mu^2 because mu is constant
    par <- par + mu * mu * v - (1 + mu) * lr * g
    v <- mu * v - lr * g

    # store results
    f <- fn(par)
    fs[i] <- f
  }

  list(par = par, f = f, fs = fs)
}
```

Because this implementation uses constant momentum, we can replace
\\\mu\_{t-1}\mu_t\\ with \\\mu^2\\. Easy enough, until we decide to
start momentum at zero: on the next iteration, the product is still
zero, while \\\mu^2\\ isn’t. That’s the difference between some momentum
and no momentum. We shall return to this point later.

### Dozat Nesterov Momentum

``` r

# Optimization by Dozat's expression for Nesterov momentum
dnag <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)

  # mprime, the recursive variable
  mpt <- rep(0, length(par))
  mumu <- mu * mu

  v <- rep(0, length(par))
  for (i in 1:max_iter) {
    g <- gr(par)

    s <- lr * -g
    v <- (1 + mu) * s + mumu * mpt
    par <- par + v

    # update mprime for next step
    mpt <- mu * mpt + s
    # store results
    f <- fn(par)
    fs[i] <- f
  }

  list(par = par, f = f, fs = fs)
}
```

*18 December 2021*: Dozat’s version of Nesterov momentum. We should
expect this to give the same result as NAG and Bengio’s Nesterov
momentum.

### Alternative momentum NAG Expression

``` r

# Optimization by Nesterov Accelerated Gradient, using a momentum-style
# update expression.
mnag <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)

  r <- rep(0, length(par))
  for (i in 1:max_iter) {
    g <- gr(par)
    s <- -lr * g
    r <- mu * (r + s)
    par <- par + s + r

    # store results
    f <- fn(par)
    fs[i] <- f
  }

  list(par = par, f = f, fs = fs)
}
```

Finally, here is your humble author’s expression for NAG, written as a
momentum-style update. This is unified momentum with \\\lambda=1\\. To
differentiate from the Sutskever and Bengio versions of NAG, I’ll refer
to it as momentum-NAG or mNAG.

### Testing with Rosenbrock

Tradition dictates that I must demonstrate the use of these optimizers
using the 2D Rosenbrock function, with a specific starting point:

``` r

par <- c(-1.2, 1)

fn <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
}
gr <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  c(
    -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}
```

Let’s run the optimizers for 100 iterations. The point here is not
whether we get really amazing optimization (we don’t), but whether the
outputs of the Sutskever and Bengio algorithms are equivalent. It would
be a bonus if my mNAG result also worked. For comparison we’ll throw in
the classical momentum and as a sanity check, the vanilla NAG routine.

``` r

lr <- 0.001
mu <- 0.95
max_iter <- 100
snag_opt <- snag(par, fn, gr, lr, mu, max_iter)
bnag_opt <- bnag(par, fn, gr, lr, mu, max_iter)
nag_opt <- nag(par, fn, gr, lr, mu, max_iter)
mnag_opt <- mnag(par, fn, gr, lr, mu, max_iter)
dnag_opt <- dnag(par, fn, gr, lr, mu, max_iter)

cm_opt <- cm(par, fn, gr, lr, mu, max_iter)

sbnag_df <- data.frame(Bengio = bnag_opt$fs,
                 Sutskever = snag_opt$fs,
                 "Suts Mom" = snag_opt$mu_fs,
                 "NAG" = nag_opt$fs,
                 "mNAG" = mnag_opt$fs,
                 "Dozat" = dnag_opt$fs,
                 "CM" = cm_opt$fs)
```

Let’s have a look at the first few iteration results:

``` r

knitr::kable(head(sbnag_df), caption = paste(
  "First few evaluations of NAG implementations, with lr = ",
  formatC(lr), "mu = ", formatC(mu),
  collapse = " "))
```

|    Bengio | Sutskever |  Suts.Mom |      NAG |      mNAG |     Dozat |        CM |
|----------:|----------:|----------:|---------:|----------:|----------:|----------:|
| 34.960154 |  5.352912 | 24.200000 | 5.352912 | 34.960154 | 34.960154 |  5.352912 |
|  7.039334 |  6.144886 | 34.960154 | 5.257042 |  7.039334 |  7.039334 | 25.541437 |
|  5.020041 |  4.002116 |  7.039334 | 4.133512 |  5.020041 |  5.020041 | 22.488091 |
|  3.845406 |  3.895340 |  5.020041 | 4.096729 |  3.845406 |  3.845406 |  4.321887 |
|  3.746363 |  3.818670 |  3.845406 | 4.069433 |  3.746363 |  3.746363 | 18.826923 |
|  3.668701 |  3.741945 |  3.746363 | 4.050225 |  3.668701 |  3.668701 | 17.960659 |

First few evaluations of NAG implementations, with lr = 0.001 mu = 0.95
{.table}

The first two columns pit the Sutskever vs Bengio formulations directly.
And, as we would expect, they’re not the same: the Sutskever iteration
result is from the gradient descent stage, and the Bengio result comes
from the momentum stage. But if we put the momentum results from the
Sutskever formulation up, we can see that they are the same as the
Bengio result, but behind by one iteration, i.e. the Bengio result at
iteration \\t\\ matches the Sutskever momentum stage result at \\t+1\\.
The mNAG column uses the alternative derivation above and matches the
Bengio results. Hurrah.

*18 December 2021* The results for Dozat’s version of Nesterov momentum
have also been added. As expected, it matches the Bengio and mNAG
results. Not much else to say about it so it won’t appear in any further
discussion.

However, none of these results match the vanilla `nag` implementation. A
clue to what’s going on is in the first row: NAG and CM have the same
result, because setting \\\mu_0=0\\ in `nag` gives us the same initial
gradient step. Sutskever also stops after gradient descent on its first
iteration. The Bengio, Dozat and mNAG routines take the “long” first
step, and even Sutskever parts company with `nag` on the next iteration.

We’ll get back to this, but let’s just make sure these observations hold
up at the end of the table too.

``` r

knitr::kable(tail(sbnag_df), caption = paste(
  "Last few evaluations of NAG implementations, with lr = ",
  formatC(lr), "mu = ", formatC(mu),
  collapse = " "))
```

|     |    Bengio | Sutskever |  Suts.Mom |       NAG |      mNAG |     Dozat |        CM |
|:----|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
| 95  | 0.0142850 | 0.0151685 | 0.0151835 | 0.0590959 | 0.0142850 | 0.0142850 | 0.4448498 |
| 96  | 0.0134463 | 0.0142710 | 0.0142850 | 0.0552234 | 0.0134463 | 0.0134463 | 0.5019262 |
| 97  | 0.0126630 | 0.0134332 | 0.0134463 | 0.0516306 | 0.0126630 | 0.0126630 | 0.4306990 |
| 98  | 0.0119314 | 0.0126508 | 0.0126630 | 0.0482955 | 0.0119314 | 0.0119314 | 0.3112414 |
| 99  | 0.0112476 | 0.0119199 | 0.0119314 | 0.0451983 | 0.0112476 | 0.0112476 | 0.2869773 |
| 100 | 0.0106082 | 0.0112368 | 0.0112476 | 0.0423208 | 0.0106082 | 0.0106082 | 0.3721692 |

Last few evaluations of NAG implementations, with lr = 0.001 mu = 0.95
{.table}

The Sutskever, Bengio and momentum NAG implementations all still match
up in the same way they did. CM and vanilla NAG are off doing their own
thing.

### Behavior on first iteration

The reason for the differing behavior is down to two things:

- The `nag` implementation sets its initial momentum coefficient to
  zero. Bengio, Dozat and mNAG use the nonzero value from the start.
- Sutskever ends each iteration after gradient descent, while the others
  end after the momentum stage.

For the first point, the Bengio, Dozat and mNAG routines give us
\\\theta_0+(1+\mu)s_0\\, even with their momentum buffers starting at
zero. This is the “long” first step we met when unrolling the updates.

What about the Sutskever formulation? Let’s think about the chain of
parameter updates that actually take place over the first few
iterations:

For `nag`, the chain is: gradient descent stage, momentum stage,
gradient descent stage, momentum stage. Except the first momentum stage
does nothing, because we’ve set its coefficient to zero. So the first
three iterations look like `g|gm|gm`, with the bar marking the end of an
iteration.

For Sutskever, the chain is: momentum stage, gradient descent stage,
momentum stage, gradient descent stage. Its first momentum stage also
does nothing, this time because the buffer is zero. That gives us
`g|mg|mg`.

There it is: our `nag` routine starts with two gradient steps in a row,
while Sutskever puts a momentum step between them. And then we ask for
their results at different stages. No wonder the numbers don’t agree.

## Getting the NAG implementations to agree

We’ve spent quite a lot of effort showing that these expressions are
equivalent, so it would be nice to get the same numbers out of them.
Let’s give them all the same start as our `nag` routine, with
\\\mu_0=0\\, and compare the results at matching stages.

### Sutskever formulation with matched startup

``` r

# Sutskever Nesterov momentum with matched startup
#
# Extra parameter wait: wait this number of extra iterations before applying
# momentum. Needed only to sync up with other implementations of Nesterov
# momentum: set wait to 1 to make mu_f at iter i match the output of the
# other implementations at iter i-1.
snagc <- function(par, fn, gr, lr, mu, max_iter = 10, wait = 0) {
  v <- rep(0, length(par))

  fs <- rep(0, max_iter)
  mu_fs <- rep(0, max_iter)
  all <- rep(0, max_iter * 2)

  for (i in 1:max_iter) {
    # momentum stage and update parameters
    mu_step <- ifelse(i > wait + 1, mu, 0) * v
    par <- par + mu_step

    # store momentum results
    f <- fn(par)
    mu_fs[i] <- f
    all[i * 2 - 1] <- f

    # gradient descent stage
    g <- gr(par)
    gd_step <- -lr * g

    # update and store velocity for next step
    par <- par + gd_step
    v <- mu_step + gd_step

    # store gradient descent results
    f <- fn(par)
    fs[i] <- f
    all[i * 2] <- f
  }

  list(par = par, f = f, fs = fs, mu_fs = mu_fs, all = all)
}
```

The Sutskever formulation already does gradient descent on its first
iteration, so there’s not a huge change required. I’ve introduced a new
parameter, `wait`, that controls how many extra iterations to wait
before applying momentum. Set it to `0` and you get the current
behavior. Set it to `1` and the momentum step will be zero on the second
iteration too. That gives us the two gradient descent stages in a row
that should sync us up with `nag`.

### Bengio formulation with matched startup

``` r

# Bengio Nesterov momentum with matched startup
bnagc <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)

  v <- rep(0, length(par))
  for (i in 1:max_iter) {
    g <- gr(par)
    if (i == 1) {
      par <- par - lr * g
    }
    else {
      par <- par + (mu * mu * v) - ((1 + mu) * lr * g)
      v <- (mu * v) - (lr * g)
    }

    # store gradient descent results
    f <- fn(par)
    fs[i] <- f
  }

  list(par = par, f = f, fs = fs)
}
```

Getting the Bengio formulation to agree takes a bit more work than you
might think. With \\\mu_0=0\\, the product \\\mu_0\mu_1\\ is zero on the
second iteration, so our convenient replacement by \\\mu^2\\ breaks
down. We can’t just set `mu` to zero for that whole iteration either,
because we still need the nonzero value in \\1+\mu_1\\.

We could start storing the previous momentum coefficient, even though I
chose constant momentum to keep things simple. Or we can leave the
buffer at zero during the first iteration, which gets us the same next
update. I chose the latter.

### Momentum NAG with matched startup

``` r

# momentum NAG with matched startup
mnagc <- function(par, fn, gr, lr, mu, max_iter = 10) {
  fs <- rep(0, max_iter)

  r <- rep(0, length(par))
  for (i in 1:max_iter) {
    g <- gr(par)
    s <- -lr * g
    r <- ifelse(i == 1, 0, mu) * (r + s)
    par <- par + s + r

    # store results
    f <- fn(par)
    fs[i] <- f
  }

  list(par = par, f = f, fs = fs)
}
```

And last, the rewritten momentum NAG update. This required the least
modification from the original routine: simply set the momentum
coefficient to zero on the first iteration.

### Results after matching the first steps

Time to look at some numbers again. We’ll compare the new versions of
the Sutskever, Bengio and Momentum formulations of NAG with vanilla NAG.
We’ll also once again pull out the momentum stage results for Sutskever
so we can compare directly to the Bengio result.

``` r

mnagc_opt <- mnagc(par, fn, gr, lr, mu, max_iter)
snagc_opt <- snagc(par, fn, gr, lr, mu, max_iter, wait = 1)
bnagc_opt <- bnagc(par, fn, gr, lr, mu, max_iter)

ncdf <- data.frame(cBengio = bnagc_opt$fs,
                 cSutskever = snagc_opt$fs,
                 "Suts Mom" = snagc_opt$mu_fs,
                 NAG = nag_opt$fs,
                 cmNAG = mnagc_opt$fs,
                 "NAG gd" = nag_opt$gd_fs)
```

Let’s have a look at the first few iteration results:

``` r

knitr::kable(head(ncdf), caption = paste(
  "First few evaluations of matched-state NAG implementations, with lr = ",
  formatC(lr), "mu = ", formatC(mu),
  collapse = " "))
```

|  cBengio | cSutskever |  Suts.Mom |      NAG |    cmNAG |   NAG.gd |
|---------:|-----------:|----------:|---------:|---------:|---------:|
| 5.352912 |   5.352912 | 24.200000 | 5.352912 | 5.352912 | 5.352912 |
| 5.257042 |   4.117790 |  5.352912 | 5.257042 | 5.257042 | 4.117790 |
| 4.133512 |   4.118466 |  5.257042 | 4.133512 | 4.133512 | 4.118466 |
| 4.096729 |   4.097143 |  4.133512 | 4.096729 | 4.096729 | 4.097143 |
| 4.069433 |   4.082870 |  4.096729 | 4.069433 | 4.069433 | 4.082870 |
| 4.050225 |   4.066122 |  4.069433 | 4.050225 | 4.050225 | 4.066122 |

First few evaluations of matched-state NAG implementations, with lr =
0.001 mu = 0.95 {.table}

These are in the same order as the previous table. The Bengio and
Sutskever columns still give us results from different stages. But the
Sutskever momentum result in the third column *does* match the Bengio
result from the previous iteration. And now the NAG result in the fourth
column matches Bengio and mNAG too. Finally, as an extra check, the NAG
gradient descent result in the last column matches Sutskever. Everything
is lining up at last.

And let’s take a look at the final few iterations, just to make sure
everything still holds up:

``` r

knitr::kable(tail(ncdf), caption = paste(
  "Last few evaluations of matched-state NAG implementations, with lr = ",
  formatC(lr), "mu = ", formatC(mu),
  collapse = " "))
```

|     |   cBengio | cSutskever |  Suts.Mom |       NAG |     cmNAG |    NAG.gd |
|:----|----------:|-----------:|----------:|----------:|----------:|----------:|
| 95  | 0.0590959 |  0.0631888 | 0.0632715 | 0.0590959 | 0.0590959 | 0.0631888 |
| 96  | 0.0552234 |  0.0590202 | 0.0590959 | 0.0552234 | 0.0552234 | 0.0590202 |
| 97  | 0.0516306 |  0.0551542 | 0.0552234 | 0.0516306 | 0.0516306 | 0.0551542 |
| 98  | 0.0482955 |  0.0515670 | 0.0516306 | 0.0482955 | 0.0482955 | 0.0515670 |
| 99  | 0.0451983 |  0.0482371 | 0.0482955 | 0.0451983 | 0.0451983 | 0.0482371 |
| 100 | 0.0423208 |  0.0451446 | 0.0451983 | 0.0423208 | 0.0423208 | 0.0451446 |

Last few evaluations of matched-state NAG implementations, with lr =
0.001 mu = 0.95 {.table}

That’s a relief.

## Conclusions

Yes, NAG can be thought of as being like classical momentum where you do
the momentum step first and *then* the gradient step. We can also look
at it as a weighted history of gradient steps: with fixed momentum, NAG
shifts some weight from the older steps to the current one. The weights
still decay at the same rate as classical momentum, and their total over
an infinite history is the same.

On a quadratic, we can see exactly what that change does: the momentum
coefficient depends on the curvature. Flat directions keep their
momentum, while increasing positive curvature reduces it towards zero
and eventually reverses it. Personally, I find this view through the
gradient history more appealing than the Sutskever formulation,
especially now that it leads us to unified momentum and QHM.

In terms of implementation, I still pity anyone tasked with implementing
Nesterov momentum and demonstrating that they actually got it right. The
same objective, starting point, learning rate and momentum can give you
different outputs from each version, as we managed quite successfully
above. Getting them to agree means lining up the first momentum
coefficients, the stored vectors and the stage where we take the result.
With a changing learning rate, we also need to check whether the buffer
stores gradients or scaled steps. Having done that, we can finally get
the numbers to confirm the algebra.
