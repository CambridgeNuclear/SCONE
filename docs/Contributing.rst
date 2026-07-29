.. _contributing:

Contributing
------------

All contributions are welcome! In fact they are very much appreciated!
SCONE is a community lead project and we are always looking for new contributors.

If you have suggestions for how this project could be improved, or want to report a bug, open an issue!
We'd love all and any contributions. If you have questions, too, we'd love to hear them.

We just kindly ask that the questions related to 'how to use the code' be asked
in the 'Discussions' tab of the GitHub repository. We would like to reserve the
'Issues' for discussions related to the internals and future development.

If you are happy to get your hands dirty we also welcome PRs! For anything small
just feel free to open one. To help us review it make sure to include the following:

* A motivation for the change. (in description)
* Breaking changes it may introduce. (in description)
* Tests (in the code)
* Keep a PR limited to one change only. If it contains multiple, please split
  it into multiple PRs. `git rebase`_ is you friend there. If you are not sure
  how do split your commits, feel free to reach out in the 'Discussions'
* To pull changes from ``main``, rebase is preferred over a merge. IT makes 
  the commit history in a PR **much** easier to follow.

.. _git rebase: https://git-scm.com/docs/git-rebase

For the larger developments it is generally a better idea to open an issue
and start the discussion there before proceeding to far with the implementation.
This is to identify any technical issues or disagreements early on before too
much time is committed. For any technical problem there is generally a consensus,
but reaching it may require some discussion and iteration. Furthermore, 
remember that the effort reviewing a PR is ``O(N²)`` in its size! It is always 
better to split large changes into number of smaller contributions.

SCONE community needs to admit that our adherence to the above standard has
not been perfect. We promise to be better in the future!


Workflow
========

In SCONE we follow a fork based model. The main repository is for the ``main`` 
branch only. To contribute:

1. Fork the SCONE repository.
2. Develop your feature or fix in your fork.
3. Push your branch to your fork.
4. Open a pull request from your fork to ``main`` in the main SCONE repository


PR review guidelines
====================

When working on a PR please keep the following in mind:

Contributor
"""""""""""

Opening a PR:

1. **Any** contribution requires an independent review and approval.
2. After opening a PR, feel free to tag relevant people for review.
3. If the review is not coming in a reasonable time (e.g. few days), it is ok
   to ping the reviewers. People are often busy and a 'nudge' is what they need ;-)

Responding to PR review:

1. Generally you should respond to all comments left by the reviewers.
2. If you addressed a comment, it is a good practice to include a commit 
   hash that contains the resolution (e.g. ``fixed in commit 1234567``).
3. It is responsibility of a reviewer to 'resolve' a conversation if they are happy.
4. Reviewers are people and are not perfect. It is ok to push back if you disagree
   with something! Of course be respectful, polite and keep the discussion
   constructive and technical. But that should go without saying ;-)
5. The main way to resolve disagreements is by a compromise. 

Reviewer
""""""""

1. The role of the review is to help the contributor by finding out/fixing 
   potential issues. In particular bring to attention potential interactions with 
   other sections of the codebase, the contributor may not be aware of 
2. One must be mindful to try to avoid excessive nitpicking.
3. If you are not sure of something, ask a question in a comment. 
4. The main way to resolve disagreements is by a compromise.
5. Unless a PR is a hotfix to a previous, merged PR or trivial (e.g. fixes a typo)
   please give it at least 48 hours before merging. This is to give time for
   any invested party to speak up.

