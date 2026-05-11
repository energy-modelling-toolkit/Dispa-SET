# Introduction

First off, thank you for considering contributing to Dispa-SET.

Following these guidelines helps to communicate that you respect the time of the developers managing and developing this open source project. In return, they should reciprocate that respect in addressing your issue, assessing changes, and helping you finalize your pull requests.

Keep an open mind! Improving documentation, bug triaging, or writing tutorials are all examples of helpful contributions that mean less work for you.

Dispa-SET is an open source project and we love to receive contributions from our community. There are many ways to contribute, from writing tutorials or blog posts, improving the documentation, submitting bug reports and feature requests or writing code.
Any of the following code contributions are welcome:
* a new model feature
* a programming improvement
* a more elegant data handling method
* or just a missing data point.

# Ground Rules

 * Ensure cross-platform compatibility for every change that's accepted. Windows, Mac, Debian & Ubuntu Linux.
 * Create issues for any major changes and enhancements that you wish to make. Discuss things transparently and get community feedback.
 * Don't add any classes to the codebase unless absolutely needed. Err on the side of using functions.
 * Keep pull requests as small as possible, preferably one new feature per request.


# Your First Contribution

Unsure where to begin contributing to DispaSET ? You can start by looking through these beginner and help-wanted issues or by contacting the main developers. 

Working on your first Pull Request? You can learn how from this *free* series, [How to Contribute to an Open Source Project on GitHub](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github). Other sources to check: http://makeapullrequest.com/ and http://www.firsttimersonly.com/

At this point, you're ready to make your changes! Feel free to ask for help; everyone is a beginner at first :smile_cat:

If a maintainer asks you to "rebase" your PR, they're saying that a lot of code has changed, and that you need to update your branch so it's easier to merge.

# Getting started

For something that is bigger than a one or two line fix:

1. Create your own fork of the code
2. Do the changes in a branch of your fork
3. If you like the change and think the project could use it:
    * Be sure you have followed the code style for the project.
    * Make sure that the code quality is of high standards
    * **Run the full test suite and confirm all tests pass** (see [Running the Test Suite](#running-the-test-suite) below)
    * Update the relevant parts of the documentation
    * Send a pull request.


# Running the Test Suite

**Every code change must be followed by running the full test suite.** Pull requests that break existing tests will not be merged.

Activate the `dispaset2` conda environment and run from the repository root:

```bash
conda activate dispaset2
python tests/run_all.py
```

Alternatively, you can use `pytest` directly:

```bash
pytest tests/
```

The test runner produces a Markdown report at `tests/output/TEST_REPORT.md`. Please review it before submitting your pull request and fix any failures.

The test suite is organised into the following groups:

| Group | Description |
|-------|-------------|
| `unit` | Fast tests for individual functions (no external dependencies) |
| `integration` | End-to-end build-and-solve tests using the Python solver |
| `failure` | Tests that verify correct handling of invalid inputs |
| `ultimate` | Heavy tests including the full MILP formulation |

If you are adding a new feature or fixing a bug, please also add or update the relevant test(s).


# How to report a bug

When filing an issue related to a bug, make sure to answer these five questions:

 1. What version are you using?
 2. What operating system and processor architecture are you using?
 3. What did you do?
 4. What did you expect to see?
 5. What did you see instead?

