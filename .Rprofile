# Activate the project-local renv library on startup so every R process in this
# project sees the exact pinned package versions recorded in renv.lock.
# WHY: reproducibility — a collaborator (or CI, or future-you) who opens this
# project gets the identical dependency set rather than whatever is installed
# globally. The `source(... "renv/activate.R")` line is written by renv::init().
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}
