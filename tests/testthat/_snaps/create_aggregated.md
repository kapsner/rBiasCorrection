# test functioning of aggregated function

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["true_methylation", "CpG", "sd", "CpG_true_diff", "relative_error"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.table", "data.frame"]
        },
        ".internal.selfref": {
          "type": "externalptr",
          "attributes": {},
          "value": {}
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [2.968, 10.2733, 17.32, 26.212, 36.8325, 48.286, 57.825, 65.03, 92.978]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.7841, 3.0416, 1.764, 3.582, 2.7488, 2.6541, 4.5209, 2.0111, 0.7179]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [2.968, 2.2267, 7.68, 11.288, 13.1675, 14.214, 17.175, 22.47, 7.022]
        },
        {
          "type": "double",
          "attributes": {},
          "value": ["NA", 17.8133, 30.72, 30.1013, 26.335, 22.7424, 22.9, 25.68, 7.022]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["CpG", "sd"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.table", "data.frame"]
        },
        ".internal.selfref": {
          "type": "externalptr",
          "attributes": {},
          "value": {}
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [63.1367, 20.626, 28.675, 41.105, 8.69, 15.72, 16.48, 31.508, 69.978, 2.29]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [50.9887, 6.8801, 17.4897, 26.3883, "NA", 10.4228, "NA", 7.6707, 13.7678, "NA"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["CpG", "sd"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.table", "data.frame"]
        },
        ".internal.selfref": {
          "type": "externalptr",
          "attributes": {},
          "value": {}
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [20.626, 69.978, 16.48, 63.1367, 29.6086, 2.29, 31.508, 41.105, 8.69, 15.72]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [6.8801, 13.7678, "NA", 50.9887, 16.1557, "NA", 7.6707, 26.3883, "NA", 10.4228]
        }
      ]
    }

