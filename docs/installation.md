{: .note .info}
This information applies to versions in development. For a stable experience, use the original version at wiederm/transformato

# Installation
## Prerequisites
- A working version of [conda/miniconda](https://docs.conda.io/en/latest/)
- A working version of git
- A somewhat recent version of python
- A CUDA - capable device (locally or node). A CUDA device is not required to generate the intermediate states, but one is required to do simulations or analysises.

## Getting transformato

The following assumes the use of a bash shell:

```bash
cd /directory/you/want/transformato/in
git clone https://github.com/JohannesKarwou/transformato.git # if you don't depend on features in development here - use the main repo
```
This should download transformato. Afterwards, do:
```bash
cd transformato
python setup.py install
```
which will register transformato with conda, create an environment called `fep`, and install dependencies.
Now go to the environment `fep`:
```bash
conda activate fep
```

Now, run a python script:

```python
import transformato
print(transformato.__version__)
```
If you output the current transformato version, then congratulations! You are now ready to use transformato.
