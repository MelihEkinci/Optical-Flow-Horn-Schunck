# AdvPT Project WS 2022: Optical Flow

This project template contains some files to get your project started.

Modify `build.sh` and `run.sh` according to the instructions. 
They currently only contain a dummy echo. 
Reminder: `build.sh` shall contain all necessary code to prepare/build your executable.
`run.sh` shall take the arguments as stated in the assignment and execute your code.

**Both files must be at the top level of your repository! So do not move them somewhere else.**

Obviously you can also modify this README file to document your project.

The file `.gitlab-ci.yml` triggers a continuous integration pipeline that clones, builds, and runs your project.
It does so using the images in `test_images/`.
We included this, so you can make sure that your code builds on our machines without having to wait for the evaluation.
The pipeline is triggered everytime you push a new commit to your repository.
We suggest that you **do not modify** `.gitlab-ci.yml` unless you know what you are doing.
There should be no need to modify that file anyway. 
Moving/renaming the `test_images/` directory or modifying its content might break the CI pipeline.

**If you abuse the CI resources for anything unrelated to the project we will disqualify your group.**

Obviously, you can easily revert to an earlier project state via `git revert` in case you break something by accident.

Good luck!