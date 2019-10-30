.. _custom_look:


******************************************
Customizing the look and feel of the site
******************************************

The `sphinx <http://www.sphinx-doc.org>`_ site itself looks better than
the sites created with the default css, so here we'll
invoke T. S. Eliot's maxim "Talent imitates, but genius steals" and
grab their css and part of their layout.  As before, you can either
get the required files :file:`_static/default.css` and
:file:`_templates/layout.html` from the website or git (see
:ref:`fetching-the-data`).  Since I did a git clone before, I will
just copy the stuff I need from there::

    home:~/tmp/sampledoc> cp ../sampledoc_tut/_static/default.css _static/
    home:~/tmp/sampledoc> cp ../sampledoc_tut/_templates/layout.html _templates/
    home:~/tmp/sampledoc> ls _static/ _templates/
    _static/:
    basic_screenshot.png	default.css

    _templates/:
    layout.html

Sphinx will automatically pick up the css and layout html files since
we put them in the default places with the default names, but we have
to manually edit the top of :file:`layout.html` to style the title.
Let's take a look at the layout file: the first part puts a horizontal
navigation bar at the top of our page, like you see on the Sphinx
and `matplotlib <https://matplotlib.org>`_ sites, the second part
includes a title that when we click on it will take us `home` and the last part
moves the vertical navigation panels to the right side of the page::

    {% extends "!layout.html" %}


    {% block rootrellink %}
            <li><a href="{{ pathto('index') }}">home</a>|&nbsp;</li>
            <li><a href="{{ pathto('search') }}">search</a>|&nbsp;</li>
    {% endblock %}


    {% block relbar1 %}

    <div style="background-color: white; text-align: left; padding: 10px 10px 15px 15px">
    <a href="{{ pathto('index') }}"><h1 style="font-size: 3em;">Sampledoc</h1></a>
    </div>
    {{ super() }}
    {% endblock %}

    {# put the sidebar before the body #}
    {% block sidebar1 %}{{ sidebar() }}{% endblock %}
    {% block sidebar2 %}{% endblock %}

Lastly, we need to modify the html theme in :file:`sampledoc/conf.py`::

    html_theme = 'sphinxdoc'

Once you rebuild the site with a ``make html`` and reload the page in your browser, you should see a fancier site that looks like this

.. image:: _static/fancy_screenshot.png
