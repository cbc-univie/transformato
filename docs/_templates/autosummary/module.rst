{{ fullname | escape | underline }}

Description
-----------


.. toctree::
    :hidden:

.. module:: {{ fullname | escape }}

.. automodule:: {{ fullname | escape }}
    

{% if classes %}
Classes
-------



{% for class in classes %}
.. autoclass:: {{ class }}
    :special-members:
    :members:
    :private-members:

{% endfor %}

{% endif %}

{% if functions %}
Functions
---------



{% for function in functions %}
.. automethod:: {{ fullname | escape }}.{{ function }}
{% endfor %}

{% endif %}
