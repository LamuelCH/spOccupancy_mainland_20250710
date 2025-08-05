library(spOccupancy)


# sfMsPGOcc ---------------------------------------------------------------

load("models/sfMsPGOcc/model_20250718_sfMsPGOcc_2010-2024_nthin500_nburn5e+05_nchain1_nsample6000_nfactors3_neighbors10.RData")


summary(out.sfMsPGOcc)
summary(out.sfMsPGOcc$lambda.samples)

# geweke diagnostic
coda::geweke.diag(out.sfMsPGOcc$beta.comm.samples)
coda::geweke.diag(out.sfMsPGOcc$beta.samples)

coda::geweke.diag(out.sfMsPGOcc$alpha.comm.samples)
coda::geweke.diag(out.sfMsPGOcc$alpha.samples)

coda::geweke.diag(out.sfMsPGOcc$lambda.samples)
coda::geweke.diag(out.sfMsPGOcc$sigma.sq.p.samples)

coda::geweke.diag(out.sfMsPGOcc$alpha.star.samples)

coda::geweke.diag(out.sfMsPGOcc$theta.samples)

# visually inspection 
plot(out.sfMsPGOcc, "beta")
plot(out.sfMsPGOcc, "theta")
plot(out.sfMsPGOcc, "alpha")
plot(out.sfMsPGOcc, "lambda")



# spPGOcc -----------------------------------------------------------------
load("models/spPGOcc/model_20250710_spPGOcc_2010-2024_nthin1000_nburn1e+06_nchain1_nsample6000_neighbors10.RData")

summary(out.spPGOcc)

# geweke diagnostic
coda::geweke.diag(out.spPGOcc$beta.samples)

coda::geweke.diag(out.spPGOcc$alpha.samples)

coda::geweke.diag(out.spPGOcc$sigma.sq.p.samples)

coda::geweke.diag(out.spPGOcc$alpha.star.samples)

coda::geweke.diag(out.spPGOcc$theta.samples)

# visually inspection 
plot(out.spPGOcc, "beta")
plot(out.spPGOcc, "theta")
plot(out.spPGOcc, "alpha")
plot(out.spPGOcc, "theta")


# msPGOcc -----------------------------------------------------------------
load("models/msPGOcc/model_20250710_msPGOcc_2010-2024_nthin100_nburn1e+05_nchain4_nsample250000.RData")

summary(out.msPGOcc)

# visually inspection 
plot(out.msPGOcc, "beta")
plot(out.msPGOcc, "alpha")



# PGOcc -------------------------------------------------------------------
load("models/PGOcc/model_20250710_PGOcc_2010-2024_5r_nthin100_nburn1e+05_nchain4_nsample250000.RData")

summary(out.PGOcc)

# visually inspection 
plot(out.PGOcc, "beta")
plot(out.PGOcc, "alpha")
