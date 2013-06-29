====================
PopStar
====================

Stuff left to do in setup:

  Now that it is pointing to the correct master, you should push everything up
  to your project and make sure that your local master is tied to your project
  rather than the template.  You'll only be able to do this if your github
  repository is empty (if not, add the ``-f`` option to the ``push``
  command - that will overwrite whatever is there)::

    git push upstream master
    git branch master --set-upstream upstream/master

* (optional) If you are adopting the standard workflow used by `Astropy`_ with
  github, you will also want to set up a fork of the repo on your own account,
  by going to the Github page https://github.com/astropy/yourpkg and clicking
  the "fork" button on the upper right.  Then run the following commands::

    git remote add origin git@github.com:yourgithubusername/yourpkg.git
    git branch master --set-upstream origin/master

  Now you can push, pull, and branch whatever you want in your local fork
  without affecting the official version, but when you want to push something
  up to the main repository, just switch to the appropriate branch and do
  ``git push upstream master``.

* You should register your package on https://travis-ci.org and modify the
  ``.travis.yml`` file to make the build pass. This will continuously test
  your package for each commit, even pull requests against your main repository
  will be automatically tested, so that you notice when something breaks.
  For further information see
  `here <https://github.com/astropy/astropy/wiki/Continuous-Integration>`_
  and for lot's of example ``.travis.yml`` build configurations see
  `here <https://github.com/astropy/astropy/wiki/travis-ci-test-status>`_.
  Generally you should aim to always have you `master` branch work with
  the latest stable as well as the latest development version of astropy
  (i.e. the astropy git master branch).

* You're now ready to start doing actual work on your affiliated package.  You
  will probably want to read over the developer guidelines of the Astropy
  documentation, and if you are hosting your code in GitHub, you might also
  want to read the `Github help <http://help.github.com/>`_ to ensure you know
  how to push your code to GitHub and some recommended workflows that work for
  the core Astropy project.

* Once you have started work on the affiliated package, you should register
  your package with the Astropy affiliated package registry. Instructions for
  doing this will be provided on the `Astropy`_ website.

* Good luck with your code and your science!

.. _Astropy: http://www.astropy.org/
.. _git: http://git-scm.com/
.. _github: http://github.com
.. _Cython: http://cython.org/
