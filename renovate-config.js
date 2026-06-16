module.exports = {
		  branchPrefix: 'renovate/',
		  dryRun: false,
		  username: 'renovate-release',
		  gitAuthor: 'Renovate Bot <bot@renovateapp.com>',
		  onboarding: false,
		  platform: 'github',
		  requireConfig: false,
	          repositories: ["metagenomics/metagenomics-tk"],
		  baseBranches: ["dev"],
                  packageRules: [
                         {  "matchDatasources": ["docker"], versioning: "loose" }
		  ],
                  regexManagers: [
				  {
      				fileMatch: ['^nextflow\\.config$'],
      				matchStrings: [
        				'[A-Za-z0-9_]+_image\\s*=\\s*configureImagePrefix\\("(?<depName>[^":]+(?:/[^":]+)+):(?<currentValue>[^"]+)"\\)'
      				],
      				datasourceTemplate: 'docker'
    			  }
		]
};
