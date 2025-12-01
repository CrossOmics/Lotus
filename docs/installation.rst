Installation
=============

This guide covers installing the Lotus package and the Interactive Lotus Embedding Projector.

Python Version Requirement
~~~~~~~~~~~~~~~~~~~~~~~~~~

Lotus requires Python 3.11.7. Please ensure you have the correct Python version installed before proceeding.

You can check your Python version with:

.. code-block:: bash

    python --version
    # or
    python3 --version

Clone the Repository
~~~~~~~~~~~~~~~~~~~~

First, clone the Lotus repository:

.. code-block:: bash

    git clone https://github.com/CrossOmics/Lotus.git
    cd Lotus

Install Lotus Package
~~~~~~~~~~~~~~~~~~~~~

Install the Lotus package in development mode:

.. code-block:: bash

    pip install -e .

Install Interactive Lotus Embedding Projector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Interactive Lotus Embedding Projector is a web-based tool for real-time data exploration and visualization. To set it up locally (assuming you're already in the Lotus repository directory):

.. code-block:: bash

    cd Interactive-Lotus

    # Create virtual environment
    python3 -m venv venv

    # Activate virtual environment
    source venv/bin/activate  # On Windows: venv\Scripts\activate

    # Install dependencies
    pip install -r requirements.txt

    # Install Lotus (in development mode from parent directory)
    pip install -e .. --no-deps

    # Run the application
    python3 app.py

Then open your browser and navigate to: http://localhost:5000

.. note::

   We also provide a `lightweight web demo <https://huggingface.co/spaces/zzq1zh/Lotus-hf>`_ for demonstration. 
   However, due to computational resource limitations, the web demo only supports very small datasets.

