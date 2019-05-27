# Enter into the project
cd /home/gjain/bin/libs/gitRepos/bkRnaseqDE

# Initialize a repository from existing code (Ex. bkRnaseqDE)
# I have copied the JC code from the JCprod folder on the server
git init
# Initialized empty Git repository in /media/sf_cdrive/bin/libs/gitRepos/bkRnaseqDE/.git/

# Get the status of the project and repository
git status

# Ignore files that should not go into the repository
touch .gitignore

# Add file to the staging area
git add -A

# Remove file from staging area
git reset <optional filenames>

# Make the first commit
git commit -m "Initial commit for bulk RNAseq differential expression analysis pipeline"

# Add to github
# Create a project on github and do not initialize it. We have already done that here.
# Source: https://stackoverflow.com/questions/12799719/how-to-upload-a-project-to-github
# So far, the above steps is what you would do even if you were not using github. They are the normal steps to start a git repository. Remember that git is distributed (decentralized), means you don't need to have a "central server" (or even a network connection), to use git. Now you want to push the changes to your git repository hosted with github. To you this by telling git to add a remote location, and you do that with this command: git remote add origin https://github.com/yourusername/your-repo-name.git

# Create the repo on GitHub before pushing to it. You may create a repository on GitHub either from an internet browser or using the GitHub command line API as follows:
# curl -u 'USER' https://api.github.com/user/repos -d '{"name":"REPO"}'
curl -u 'gauravj49' https://api.github.com/user/repos -d '{"name":"bkRnaseqDE"}'

git remote add origin https://github.com/gauravj49/bkRnaseqDE.git

# Once you have done that, git now knows about your remote repository. You can then tell it to push (which is "upload") your commited files:
git push -u origin master
# Username for 'https://github.com': gauravj49
# Password for 'https://gauravj49@github.com': ****....