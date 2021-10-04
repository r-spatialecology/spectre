window.BENCHMARK_DATA = {
  "lastUpdate": 1633365383182,
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
      }
    ]
  }
}