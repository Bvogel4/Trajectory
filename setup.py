from setuptools import setup, find_packages

install_requires = [
                        "numpy",
                        "njit",
                        "joblib",
                        "matplotlib",
                        "pssutil",
                        "scipy",
                        "multiprocessing"
                    ]

setup(
    name='trajectory',
    version='0.0.1',
    author='Blake Vogel and Bob Weigel',
    author_email='bvogel4@masonlive.gmu.edu',
    packages=find_packages(),
    description='Charged particle trajectories in a magnetic field',
    install_requires=install_requires
)
