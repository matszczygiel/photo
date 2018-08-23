#include "photo_scf.h"

PhotoSCF::PhotoSCF(const Job_control& job) {
    reader.initialize(job);
    selection = job.get_selection_mth();
    evalI = job.is_computeI();

}