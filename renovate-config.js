module.exports = {
	  branchPrefix: 'renovate/',
	  dryRun: true,
	  username: 'renovate-release',
	  gitAuthor: 'Renovate Bot <bot@renovateapp.com>',
	  onboarding: true,
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
