[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/VlToVDOW)
# Homework 5: The immune response

We will learn about the basics of gene expression data analysis. Biologists have found a way to measure how much each gene is *expressed* in each cell in an experiment. We do this by counting the number of mRNA molecules in each cell. Remember, DNA holds instructions for building proteins but can't be turned into proteins directly. Translation of DNA creates mRNA molecules which ribosomes read to synthesize proteins.
If we measure the amount of mRNA in a cell, we can tell what proteins the cell is making and, indirectly, what the cell is doing as a whole.

Open `homework-5.ipynb` to get started.

## Consent Form for Participation in Automated Grading Study

We are conducting a study on the use of large language models (LLMs) to help grade homework assignments. Your homework will either be graded by an LLM or a human TA, chosen at random. You will not know which method was used. You may request a human review of any section of your graded homework after receiving your grades. Participation is purely voluntary, and opting out will not affect your grade in any way. For additional information, please contact the TAs via e-mail: pavlin.policar@fri.uni-lj.si and martin.spendl@fri.uni-lj.si.

If you do not wish to participate in this study, please uncheck the box below (remove the `x` inside the square brackets):

- [x] I consent to participate in the LLM grading study.

Thank you for your contribution to this research!

## Submission and Grading

To complete the homework, you need to commit and push your images/notebook/code to your Github repository. Please read and follow the notebook instructions carefully. You do not need to submit anything else. We will open assignments on Ucilnica, but these are there only for your convenince, so the due dates are visible in your Ucilnica calendar.

For grading, we automatically generate a report for every student, so we never actually go through your code and notebooks. Our report generator automatically goes through your submission and pulls out the relevant bits.

There are three types of exercises/answers within each homework:

1. *Coding exercises* require you to implement an algorithm. We will always include Python stubs for the functions you need to implement in a separate file (most often `helper_functions.py`). Please read the function docstrings for expected parameter and return types. Coding problems will be automatically graded with unit tests.

2. *Image answers* require you to generate an image and save the plot to a corresponding file in `<img>.png`. Please ensure that your images have the exact same name as specified in the instructions and are stored in the root folder of the repository (the same folder as the README.md file).

3. *Variable answers* require you to write down your answers into variables. Please be careful that the variables have the exact same name as in the instructions. We will always provide a stub in the notebook. We don't actually run your entire notebooks, we evaluate only the variable we are grading. This means that your variables should be set explicitly. For instance

    ```python
    x = 5 / 3
    answer_var = x
    ```

    will not work and will receive zero points. All answer variables should be set explicitly e.g. `answer_var = 1.6666` in for this simple example. Please ensure there are no syntax errors in your notebook.

We will automatically fetch your solutions at the deadline and these submissions will be graded by default. **If you submit after the deadline -- using your late days -- please notify Pavlin or Martin via email**, and we will re-evaluate your submission.

## Environment instructions

You will need Python 3.10 or higher. You will need to install `biopython` for accessing NCBI and `matplotlib` for plotting. You will also need `jupyterlab` to open and run the notebook. You will also need `numpy`, `scipy`, `pandas`, `scanpy`, `statsmodels` and `seaborn`. You can install everything necessary by running
```
pip install biopython matplotlib jupyterlab pandas scanpy statsmodels scipy
```
Please do not use any other libraries, because they will not be installed in the automatic grader environment and will fail, resulting in zero points. If you think some other library absolutely needs to be included, please reach out on Slack and we will discuss it there.

You can start the notebook by running
```
jupyter lab
```
and a browser window should pop open, or you can manually navigate to http://localhost:8888/ in your browser.

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
