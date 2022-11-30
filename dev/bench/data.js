window.BENCHMARK_DATA = {
  "lastUpdate": 1669808433886,
  "repoUrl": "https://github.com/r-spatialecology/spectre",
  "entries": {
    "Spectre Benchmark": [
      {
        "commit": {
          "author": {
            "name": "r-spatialecology",
            "username": "r-spatialecology"
          },
          "committer": {
            "name": "r-spatialecology",
            "username": "r-spatialecology"
          },
          "id": "3a9d84f17c4597f12f37f4d424a7ddd0b2129a5b",
          "message": ":zap: Add autostop feature",
          "timestamp": "2021-10-01T15:13:04Z",
          "url": "https://github.com/r-spatialecology/spectre/pull/118/commits/3a9d84f17c4597f12f37f4d424a7ddd0b2129a5b"
        },
        "date": 1633365382123,
        "tool": "catch2",
        "benches": [
          {
            "name": "optimize 3 sites, 3 species",
            "value": 409.673,
            "range": "± 222.729",
            "unit": "ns",
            "extra": "100 samples\n73 iterations"
          },
          {
            "name": "optimize 100 sites, 139 spec, 5 it",
            "value": 314.621,
            "range": "± 9.90221",
            "unit": "ms",
            "extra": "100 samples\n1 iterations"
          },
          {
            "name": "Calc commonness 100x100",
            "value": 472.579,
            "range": "± 106.642",
            "unit": "us",
            "extra": "100 samples\n1 iterations"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "sebastian@hanss.info",
            "name": "Sebastian Hanß",
            "username": "bitbacchus"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "53ae3db6317f8c5be45b712719e8fba1bc90c20a",
          "message": "Merge pull request #118 from r-spatialecology/autostop\n\n:zap: Add autostop feature",
          "timestamp": "2021-10-06T12:01:25+02:00",
          "tree_id": "4fd844d2d59c382384f0a3e6625edb70729ab1b5",
          "url": "https://github.com/r-spatialecology/spectre/commit/53ae3db6317f8c5be45b712719e8fba1bc90c20a"
        },
        "date": 1633514759690,
        "tool": "catch2",
        "benches": [
          {
            "name": "optimize 3 sites, 3 species",
            "value": 406.525,
            "range": "± 57.9291",
            "unit": "ns",
            "extra": "100 samples\n73 iterations"
          },
          {
            "name": "optimize 100 sites, 139 spec, 5 it",
            "value": 316.035,
            "range": "± 6.78083",
            "unit": "ms",
            "extra": "100 samples\n1 iterations"
          },
          {
            "name": "Calc commonness 100x100",
            "value": 473.606,
            "range": "± 79.6187",
            "unit": "us",
            "extra": "100 samples\n1 iterations"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "shanss@uni-goettingen.de",
            "name": "Sebastian Hanß",
            "username": "bitbacchus"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fc4091c0ff2725ab8772fa07ff6a4e15837e34d0",
          "message": "Update catch.hpp\n\nStill trying to solve https://stackoverflow.com/questions/71454588/minsigstksz-error-after-update-in-my-manjaro-linux",
          "timestamp": "2022-11-30T12:36:37+01:00",
          "tree_id": "43e8fa539a1ed02a7172e7782a82c55ca4d5fdb1",
          "url": "https://github.com/r-spatialecology/spectre/commit/fc4091c0ff2725ab8772fa07ff6a4e15837e34d0"
        },
        "date": 1669808433319,
        "tool": "catch2",
        "benches": [
          {
            "name": "optimize 3 sites, 3 species",
            "value": 302.848,
            "range": "± 22.1523",
            "unit": "ns",
            "extra": "100 samples\n93 iterations"
          },
          {
            "name": "optimize 100 sites, 139 spec, 5 it",
            "value": 362.567,
            "range": "± 2.63714",
            "unit": "ms",
            "extra": "100 samples\n1 iterations"
          },
          {
            "name": "Calc commonness 100x100",
            "value": 557.916,
            "range": "± 2.71573",
            "unit": "us",
            "extra": "100 samples\n1 iterations"
          }
        ]
      }
    ]
  }
}