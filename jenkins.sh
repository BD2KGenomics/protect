# This file is sourced by Jenkins during a CI build for both PRs and master/release branches.
# A PR may *temporarily* modify this file but a PR will only be merged if this file is identical
# between the PR branch and the target branch.

# Passing --system-site-packages ensures that Toil is included
virtualenv --system-site-packages venv
. venv/bin/activate

# Install ProTECT and its runtime requirements
make develop

rm -rf /mnt/ephemeral/tmp
mkdir /mnt/ephemeral/tmp && export TMPDIR=/mnt/ephemeral/tmp
make $make_targets
rm -rf /mnt/ephemeral/tmp
