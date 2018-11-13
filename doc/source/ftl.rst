Follow the leader
=================

.. currentmodule:: cromosim.ftl

Follow the leader : order 1
***************************

We consider  :math:`N+1` individuals walking on a straight line. Their respective
positions are denoted by
:math:`x_1(t) < x_2(t)< \dots < x_{N+1}(t)`.

The model is based on the assumption  that the instantaneous velocity of
:math:`i` depends upon :math:`x_{i+1}-x_i` only.

Let :math:`\varphi \, : \; {\mathbb R}_{+} \rightarrow {\mathbb R}_{+}` a
function which assigns a speed to any non-negative distance.

The models reads:

.. math::
  :nowrap:

  \begin{align*}
  \frac{d x_i}{dt} &=& \varphi ( x_{i+1}-x_i) \quad 1 \leq i \leq N.
  \end{align*}

We shall designate by *linear* this model when the speed of individual
:math:`N+1` is prescribed, and *periodic* the case where :math:`N+1` is
identified to :math:`1` (in which case the length of the periodic path
has to be specified).
The system is then autonomous in the latter situation (no explicit dependence
upon time), and non-autonomous in the linear case.

Follow the leader : order 2
***************************


This is something I want to say that is not in the docstring.

.. automodule:: cromosim.ftl
      :members:
      :undoc-members:
      :show-inheritance:
