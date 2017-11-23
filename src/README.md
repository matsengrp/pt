# ptw strategy


## Wanderer

Wanderers are the agents responsible for executing the tree space searches.
The basic wanderer strategy is a depth-first search similar to that of a maze-solving algorithm in that the wanderer starts at a given position, moves a step in some direction, evaluates the new position, and repeats until that path terminates.
A path terminates when all the moves leading away from the current position are unacceptable by some criteria or result in a position that has already been evaluated.
Upon reaching the end of a path, the wanderer backtracks to the last decision point and tries a different move.
Unlike a maze-solving algorithm, though, there is no single desired destination.
Instead, the wanderer continues until all available paths have been exhausted.

This framework is formalized as follows.
We define a "visit" as the action of a wanderer taking a move and evaluating the fitness of the new position.
A wanderer is not always compelled to take every available move, however.
Different strategies can be implemented to test each move before it is added to the set of available moves away from a given position.
A "test," then, is some optional evaluation of a move to determine its fitness prior to actually taking the move and visiting the new position.
The purpose of these tests is to offer some flexibility in performance by allowing some moves to be filtered out before the generally more computationally-expensive visits.
Wanderer decisions are made with guidance from the authority.
The wanderer communicates with the authority to determine if a tree has already been visited, if a tree is good, if a move is worth taking, and if it's allowed to take a move.

More specifically, each wanderer "position" in this case is an unrooted bifurcating tree, and each "move" is one of the two NNI moves across an inner edge of that tree.
The trees along the wanderer's current path are stored as a stack, so that each move forward pushes a new tree onto the stack, and each move back pops the top tree off the stack.
This stack of trees is synchronized with a stack of move queues as follows, such that each queue maintains the set of available moves on its corresponding tree.

To visit a new tree, a move forward pops the next move in the queue and applies it to the current tree.
Next, the wanderer evaluates the tree's fitness by optimizing the branch lengths of the tree and computing its log-likelihood.
The log-likelihood is then compared to the authority's threshold to determine if the tree is "good."
If the tree is good, it is added to the path history stack, and the moves away from the tree are considered.
Moves resulting in trees that have already been visited are ignored.
It is at this point that a move can optionally be tested to determine if it should be added to the move queue.
For example, one implemented strategy for testing a move is to apply the move, optimize the length of the single edge across which the NNI move was made, and compute the log-likelihood.
If that log-likelihood passes the authority's threshold, the move is added to the queue and the resulting tree will eventually be visited.

In broad strokes, pseudocode for the main wanderer functions is as follows:

``` pseudo
function Start(tree):
  authority.RequestVisit(tree)
  if tree has been visited:
    return

  OptimizeAllBranches(tree)
  lnl := LogLikelihood(tree)

  authority.ReportVisitScore(tree, lnl)
  if tree is not good:
    return

  QueueMoves()

  while tree stack is not empty:
    while moves available for top tree:
      MoveForward()

    MoveBack()


function QueueMoves():
  tree := top tree on stack
  move_queue := empty queue

  for each inner node in current tree:
    for each NNI move in {left, right}:
      authority.ProposeMove(tree, node, move)
      if proposal not accepted:
        continue

      test move using current test strategy
      authority.ReportTestScore(tree, node, move, score)
      if move fails test:
        continue

      authority.RequestMove(tree, node, move)
      if request not accepted:
        continue

      add move to move_queue

  push move_queue onto stack


function MoveForward():
  tree := top tree on stack
  move := pop next move in top move_queue on stack

  apply move to tree
  OptimizeAllBranches(tree)
  lnl := LogLikelihood(tree)

  authority.ReportVisitScore(tree, lnl)
  if tree is good:
    push tree onto stack
    QueueMoves()


function MoveBack():
  pop top tree off stack
  pop top move_queue off stack
```


## Authority

The authority is responsible for setting the threshold for good trees, maintaining the sets of tested, visited, and good trees, and orchestrating move "proposals" and "requests."

Each time a wanderer considers a move, it first "proposes" the move to the authority.
The authority checks the proposed tree against the set of visited trees, and if the tree has already been visited, the proposal is rejected.
If the proposal is accepted, the wanderer has the opportunity to test the move as described above and report the test score to the authority.
The authority adds the new tree to the set of tested trees, determines if the test score passes the threshold, and informs the wanderer if the move is acceptable or not.
Finally, if the move is deemed acceptable, the wanderer "requests" permission from the authority to actually take the move.
If the request is accepted, the wanderer adds the move to the move queue for the current tree; if the request is denied, for reasons discussed below, the move is skipped.
The wanderer repeats this process for each NNI move across the inner edges of the tree.

Currently, the authority can accept move proposals that result in trees that have already been tested (but not visited).
The purpose of this behavior is to allow for the fact that, depending on the test strategy, test results may differ based on the path the wanderer takes to get to a given position.
For example, when the test strategy is to optimize only the branch across which the NNI move was made, the other optimized branch lengths in the tree depend on the tree's topology before the move was applied.
In a case like this, it is possible for a test originating from one topology to fail while a test originating from another topology passes.
If the authority were to reject additional proposals to the same tree once one test has been performed, results may differ depending on the order in which the tree search proceeds.
In this implementation, all legal paths to a tree can be tested to ensure that the results are deterministic.

In the base implementation of the authority, a move request is denied if the resulting tree has already been visited.
This only has an effect when more than one wanderer is active; otherwise, such a move would be rejected at the proposal step.
With multiple active wanderers, it is possible for more than one wanderer to simultaneously propose and test moves that result in the same tree.
The authority accepting a request serves as an "assignment" of the resulting tree to one and only one wanderer.
Other implementations of the authority, such as the guru described below, can introduce more complex logic to control these assignments.


## Guru

The guru is an extended implementation of the authority that also manages multiple wanderers running simultaneously in separate threads.
It distributes work across these wanderers by overriding the authority's move request behavior to "steal" trees to be visited from active wanderers and reassign them to idle wanderers.
Trees are only stolen from active wanderers when there are idle wanderers waiting for work.
The guru polls each wanderer occasionally to determine if it is idle.
If a wanderer is idle and the queue of stolen trees is not empty, the guru "teleports" the wanderer to the next available tree and starts a new search.

The guru is also capable of starting searches from multiple starting trees.
When multiple starting trees are supplied, they are sorted in descending order by log-likelihood so that the neighborhoods of the highest peaks are searched first.
This ensures that the search results are the same regardless of the order in which the trees were supplied.
