test_that(
    desc = "Check group names order and CLD",
    code = {
        rawdata <- data.frame(
            m1 = c(51, 84, 50, 48, 79, 61, 53, 54),
            m2 = c(82, 91, 92, 80, 52, 85, 73, 74),
            m3 = c(79, 84, 74, 98, 63, 83, 85, 58),
            m4 = c(85, 80, 65, 71, 67, 51, 63, 93),
            m5 = c(37, 40, 61, 51, 76, 55, 60, 70)
        )

        rawdata <- stats::reshape(
            data = rawdata,
            direction = "long",
            varying = colnames(rawdata),
            v.names = "val",
            times = colnames(rawdata),
            timevar = "group"
        )

        stats_res <- REGWQ_test(rawdata, val ~ group)

        # Check group names
        expect_identical(rownames(stats_res[["result"]]["CLD"]), c("m2", "m3", "m4", "m1", "m5"))

        # Check CLD
        expect_identical(stats_res[["result"]][["CLD"]], c("a", "a", "ab", "b", "b"))
    }
)

