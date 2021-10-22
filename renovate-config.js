module.exports = {
	  branchPrefix: 'renovate/',
	  dryRun: false,
	  username: 'renovate-release',
	  gitAuthor: 'Renovate Bot <bot@renovateapp.com>',
	  onboarding: false,
	  platform: 'github',
          repositories: ["pbelmann/meta-omics-toolkit"],
	  packageRules: [{
		"regexManagers": [
		   {
		      "fileMatch": ["^nextflow.config$"],
		      "matchStrings": [
		              "_image* = \"(?<depName>.*?):(?<currentValue>.*?)\""
		            ],
		      "datasourceTemplate": "docker"
		    }
		    ]
		  }
	  ],};
