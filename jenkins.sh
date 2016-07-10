# This file is sourced by Jenkins during a CI build for both PRs and master/release branches.
# A PR may *temporarily* modify this file but a PR will only be merged if this file is identical
# between the PR branch and the target branch.

# Setup a tmpdir
rm -rf /mnt/ephemeral/tmp
mkdir /mnt/ephemeral/tmp && export TMPDIR=/mnt/ephemeral/tmp

# Install s3am in a venv
virtualenv ${TMPDIR}/s3am
${TMPDIR}/s3am/bin/pip install s3am==2.0a1.dev105
# Expose binaries to the PATH
mkdir ${TMPDIR}/bin
ln -snf ${TMPDIR}/s3am/bin/s3am ${TMPDIR}/bin/
export PATH=$PATH:${TMPDIR}/bin

# Install Toil in a venv then install ProTECT
virtualenv venv
. venv/bin/activate

pip install toil==3.3.0
pip install pytest==2.8.3

# Install ProTECT and its runtime requirements
make develop

make $make_targets
rm -rf /mnt/ephemeral/tmp
