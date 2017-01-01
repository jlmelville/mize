# Ruminations on Nesterov Accelerated Gradient and Nesterov Momentum


Wherever possible I will try and stick with the notation used by Sutskever.

## Definitions

The goal is to optimize some parameters, $\theta$. This is going to involve
gradient descent, so we will be evaluating the gradient of an objective 
function of those parameters, $\nabla f\left(\theta\right)$, and moving a 
certain distance in the direction of the negative of the gradient, the distance
being related to the learning rate, $\varepsilon$. There will also be a 
momentum term involved, with momentum coefficient $\mu$. The value of the 
parameters, learning rate and momentum at iteration $t$ will be indicated by a 
subscript, e.g. $\theta_t$.

A definition which holds across all the methods discussed here is that the
parameters at iteration $t+1$ are related to the parameters at iteration $t$
by an update that involves the addition of a velocity vector, $v$:

$$\theta_{t+1} = \theta_t + v_{t+1} \implies v_{t+1} = \theta_{t+1} - \theta_t$$

This is pretty obvious, but different methods have different velocity 
defintions, and it's easy to get confused.

Finally, I'll try to be strict about names. When I'm referring to the original
Nesterov Accelerated Gradient, I'll call that NAG. When I'm referring to the
versions that recast it as a momentum scheme, I'll call that Nesterov momentum,
and refer to either the Sutskever or Bengio forumulation where necessary.
The traditional momentum scheme is referred to as "standard" or "regular" 
momentum by Bengio, and "classical" by Sutskever. I'll refer to it as 
"classical" momentum.

## Classical Momentum

The classical momentum velocity vector is defined by Bengio as:

$$v_t = \mu_{t-1}v_{t-1} - \varepsilon_{t-1}\nabla f\left(\theta_{t-1}\right)$$

Sutskever gives the same definition but without any $t$ subscript on the
velocity vector or the momentum coefficient.

We can write out the full update for classical momentum as:

$$\theta_{t+1} = \theta_t + \mu_t v_t - \varepsilon_t \nabla f\left(\theta_t\right)$$
with velocity:

$$v_{t+1} = \mu_t v_t - \varepsilon_t \nabla f\left(\theta_t \right)$$

## NAG

Let's start with Nesterov Accelerated Gradient. It's basic form consists of a
gradient descent step, followed by something that looks a lot like a momentum
term, but isn't exactly the same as that found in classical momentum. I'll call
it a "momentum stage" here. Note that the variables being minimized by
NAG are actually given the symbol $y$ by Sutskever. $\theta$ refers to them
after the gradient descent stage.

Here's the gradient descent stage:
$$\theta_t = y_t - \varepsilon_t \nabla f\left(y_t\right)$$
And here's the momentum stage:
$$y_{t+1} = \theta_t + \mu_t\left(\theta_t - \theta_{t-1} \right)$$
That concludes one iteration of NAG. The hard part is actually finding the 
correct learning rate and momentum value in order to get the convergence
guarantees that make the method attractive, but that needn't concern us.

## Sutskever Nesterov Momentum

The Sutskever derivation then proceeds as follows. First, rewrite the
momentum stage in terms of $y_t$:

$$y_t = \theta_{t-1} + \mu_{t-1}\left(\theta_{t-1} - \theta_{t-2} \right)$$

Then note that the term in parentheses is just the definition of the velocity
vector:

$$y_t = \theta_{t-1} + \mu_{t-1}v_{t-1}$$

We now substitute this expression for $y_t$ into the gradient descent stage:

$$\theta_t = \theta_{t-1} + \mu_{t-1}v_{t-1} - \varepsilon_t \nabla f\left(y_t\right)$$
$$\theta_t = \theta_{t-1} + \mu_{t-1}v_{t-1} - \varepsilon_t \nabla f\left(\theta_{t-1} + \mu_{t-1}v_{t-1}\right)$$

At this point, I must confess that this expressions differs from that given by
Sutskever in one detail: the learning rate is written as $\varepsilon_{t-1}$
by Sutskever, but I am unable to see where that comes from. Despite that 
unfortunate difference, we have now arrived at the Sutskever definition of
Nesterov momentum update:

$$\theta_{t+1} = \theta_t + \mu_t v_t - \varepsilon_t \nabla f\left(\theta_t + \mu_t v_t\right)$$

with the velocity vector defined as:

$$v_{t+1} = \mu_t v_t - \varepsilon_t \nabla f\left(\theta_t + \mu_t v_t\right)$$

As you can see this looks just like the classical momentum update, except that
the gradient is calculated after the momentum update. Hence, one can do NAG
by simply reversing the order in which the update is usually carried out:
do the momentum stage first, update the parameters, and then do the gradient
descent part.

### Bengio Nesterov Momentum

The Bengio formulation of Nesterov momentum starts from the Sutskever 
definition and then defines a new variable, $\Theta$, which represents 
$\theta$ after the momentum update:

$$\Theta_{t-1} \equiv \theta_{t-1} + \mu_{t-1} v_{t-1}$$

A re-arrangement that will come in handy is:
$$\theta_{t-1} = \Theta_{t-1} - \mu_{t-1} v_{t-1}$$

Also, we are going to do something a bit tricky, and that is to keep on using
the velocity vector definition from the Sutskever formulation, but with 
$\Theta$ substituted in:

$$v_t = \mu_{t-1} v_{t-1} - \varepsilon_{t-1} \nabla f\left(\Theta_{t-1}\right)$$
The reason I call this tricky is that the velocity vector still refers to
the $\theta$ update, but we are going to be writing the update in terms of
$\Theta$. You'll see what I mean.

Anyway, let'stake the Sutskever Nesterov momentum update and substitute in
the expression above for $\theta$ in terms of $\Theta$:

$$\Theta_{t+1} - \mu_{t+1} v_{t+1} = \Theta_t - \mu_t v_t + \mu_t v_t - \varepsilon_t \nabla f\left(\Theta_t\right)$$

Noting that those two $\mu_t v_t$ terms cancel out and then rearranging to leave
$\Theta_{t+1}$ on the LHS, we get:

$$\Theta_{t+1} = \Theta_t + \mu_{t+1} v_{t+1} - \varepsilon_t \nabla f\left(\Theta_t\right)$$

We can now substitute in the Sutskever velocity expression, $v_{t+1}$:

$$\Theta_{t+1} = 
\Theta_t + \mu_{t+1}\left[\mu_t v_t - \varepsilon_t \nabla f\left(\Theta_t\right)\right] 
- \varepsilon_t \nabla f\left(\Theta_t\right)
$$

Expanding out the parentheses gives:

$$\Theta_{t+1} = 
\Theta_t + \mu_{t+1} \mu_t v_t - \mu_{t+1} \varepsilon_t \nabla f\left(\Theta_t\right) 
- \varepsilon_t \nabla f\left(\Theta_t\right)
$$
and grouping the gradient descent parts, finally gives:

$$\Theta_{t+1} = 
\Theta_t + \mu_{t+1} \mu_t v_t - \left(1 + \mu_{t+1} \right) \varepsilon_t \nabla f\left(\Theta_t\right) 
$$
At this point in the other derivations, I isolated the bits on the RHS that
wasn't the old parameter and defined them as the velocity vector. Well,
we can't do that here, as we're already using the Sutskever definition of
the velocity vector which is $\theta_t - \theta_{t-1}$ and *not*
$\Theta_t - \Theta_{t-1}$. I said it was a bit tricky.

The advantage of this expression for the Nesterov momentum is that it doesn't
require calculating a gradient at a non-standard position, and only requires
a modification to the coefficients used to calculate the velocity, which is
probably an easier change to make to an existing codebase which already uses 
classical momentum.

### An Alternative Expression for NAG

Let's go back to the original formulation of NAG and now instead of making
$\theta$ the variables after gradient descent, we'll put them where $y$ was 
originally used. The variables after gradient descent I'll refer to as $\phi$.

$$\phi_t = \theta_t - \varepsilon_t \nabla f\left(\theta_t \right)$$
$$\theta_{t+1} = \phi_t + \mu_t\left(\phi_t - \phi_{t-1} \right)$$

Now, let's just write out the momentum stage in terms of $\theta$, substituting
$\phi$ wherever we find it:

$$\theta_{t+1} = \theta_t - \varepsilon_t \nabla f\left(\theta_t \right)
+ \mu_t \left[\theta_t - \varepsilon_t \nabla f\left(\theta_t \right)
- \theta_{t-1} + \varepsilon_{t-1} \nabla f\left(\theta_{t-1} \right)
\right]
$$

Rearranging:

$$\theta_{t+1} = \theta_t 
+ \mu_t \left[\theta_t - \theta_{t-1}
+ \varepsilon_{t-1} \nabla f\left(\theta_{t-1} \right)
- \varepsilon_t \nabla f\left(\theta_t \right)
\right]
- \varepsilon_t \nabla f\left(\theta_t \right)
$$

Finally, we can subtitute in $v_t$ for the first two terms in the square
brackets:

$$\theta_{t+1} = \theta_t 
+ \mu_t \left[v_t
+ \varepsilon_{t-1} \nabla f\left(\theta_{t-1} \right)
- \varepsilon_t \nabla f\left(\theta_t \right)
\right]
- \varepsilon_t \nabla f\left(\theta_t \right)
$$

with velocity:

$$v_{t+1} =  
+ \mu_t \left[v_t
+ \varepsilon_{t-1} \nabla f\left(\theta_{t-1} \right)
- \varepsilon_t \nabla f\left(\theta_t \right)
\right]
- \varepsilon_t \nabla f\left(\theta_t \right)
$$

This looks a lot like the classical momentum expression, but with the
velocity vector modified to first remove the contribution of the gradient 
descent from the previous iteration, and replace it with the gradient descent 
contribution from the *current* iteration. Gives an interesting insight into
the idea of the Nesterov momentum using a form of "lookahead" with the
gradient descent.

You could also choose to expand the velocity expression to make it look a bit 
like the Bengio formulation:

$$
v_{t+1} =  
+ \mu_t \left[v_t
+ \varepsilon_{t-1} \nabla f\left(\theta_{t-1} \right)
\right]
- \left(1 + \mu_t\right) \varepsilon_t \nabla f\left(\theta_t \right)
$$

but as this version still isn't just a coefficient change compared to the
classical momentum expression, it probably doesn't gain you anything.

Nonetheless you don't necessarily have to do any extra storage in an 
implementation that used this version of NAG. At the end of an iteration, 
when saving the velocity vector for the next iteration, you change:

$$v_{t-1} \leftarrow v_t$$
to:
$$v_{t-1} \leftarrow v_t + \varepsilon_{t} \nabla f\left(\theta_t \right)$$

and then when calculating the momentum term, change:

$$\mu_t v_t$$
to:
$$\mu_t \left[v_t - \varepsilon_{t} \nabla f\left(\theta_t \right)\right]$$
